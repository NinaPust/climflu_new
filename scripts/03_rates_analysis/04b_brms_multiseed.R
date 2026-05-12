# scripts/03_rates_analysis/04b_brms_multiseed.R
#
# PURPOSE: Stability diagnostic for brms phylogenetic regression.
#          Runs the same model over N fixed seeds, each drawing a fresh
#          stratified subsample (year × phylogenetic cluster, N_SUB = 1000).
#          Pools posteriors across seeds using the law of total variance and
#          flags predictors where between-seed variance exceeds within-seed
#          variance (unstable = subsample sensitivity is larger than
#          estimation uncertainty).
#
# OUTPUT:
#   data/rates/brms/{gene}_brms_multiseed.csv   — raw per-seed posteriors
#   data/rates/brms/all_brms_multiseed_raw.csv  — all genes stacked
#   data/rates/brms/all_brms_multiseed.csv      — pooled summary
#
# FLAGS:
#   --gene X       run a single gene segment
#   --test         use subsampled data from data/test/
#   --chains N     number of MCMC chains (default 4)
#   --iter N       total iterations per chain (default 2000)

suppressPackageStartupMessages({
  library(here)
  library(brms)
  library(ape)
  library(dplyr)
  library(readr)
  library(tidybayes)
  library(cluster)
})

source(here("config.R"))
source(here("R", "climate_utils.R"))

# ── Settings ──────────────────────────────────────────────────────────────────
N_SUB  <- 1000
SEEDS  <- c(42, 123, 7, 999, 2024, 314, 888)

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

args   <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  as.numeric(args[idx + 1])
}

genes    <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}
use_test <- "--test" %in% args
CHAINS   <- as.integer(get_arg("--chains", 4))
ITER     <- as.integer(get_arg("--iter",   2000))

if (use_test) message("*** Running in TEST mode ***")

paths_timed <- if (use_test) PATHS$test_trees_timed else PATHS$trees_timed
paths_rates <- if (use_test) PATHS$test_rates       else PATHS$rates
paths_out   <- if (use_test) PATHS$test_rates       else {
  p <- file.path(PATHS$rates, "brms")
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
  p
}

PREDICTOR_TERMS <- c("scaled_temp", "scaled_humidity",
                     "scaled_delta_temp", "scaled_delta_humidity",
                     "scaled_sampling_density", "scaled_pop_density")

# ── Helper: compute phylogenetic cluster assignments (once per gene) ──────────
# Returns df with a phylo_cluster column added. Cophenetic distances are O(n^2)
# so this should be called once and reused across seeds.
add_phylo_clusters <- function(df, tree, n_clust = 40) {
  if (!requireNamespace("cluster", quietly = TRUE)) return(df)
  message("  Computing cophenetic distances for phylogenetic clustering...")
  coph    <- cophenetic(tree)
  tips    <- intersect(rownames(coph), df$accession)
  coph    <- as.dist(coph[tips, tips])
  k       <- min(n_clust, length(tips) - 1L)
  pam_fit <- cluster::pam(coph, k = k, diss = TRUE)
  df$phylo_cluster <- pam_fit$clustering[df$accession]
  df
}

# ── Helper: stratified subsample (clusters pre-computed) ─────────────────────
stratified_subsample <- function(df, n, seed = 42L) {
  set.seed(seed)
  if (nrow(df) <= n) return(df)

  strat_vars <- intersect(c("phylo_cluster", "year"), names(df))
  if (length(strat_vars) == 0) return(df[sample(nrow(df), n), ])

  strata <- df |>
    count(across(all_of(strat_vars)), name = "grp_n") |>
    mutate(.k = pmax(1L, as.integer(round(grp_n / nrow(df) * n))))

  sampled <- df |>
    left_join(strata, by = strat_vars) |>
    group_by(across(all_of(strat_vars))) |>
    group_map(function(sub, key) {
      k <- unique(sub$.k)
      sub[sample(nrow(sub), min(k, nrow(sub))), ]
    }) |>
    bind_rows()

  sampled[sample(nrow(sampled), min(n, nrow(sampled))), ]
}

# ── Helper: posterior summary for fixed effects ───────────────────────────────
posterior_summary_brms <- function(fit, predictors) {
  draws <- as_draws_df(fit)

  lapply(c("Intercept", predictors), function(term) {
    col <- paste0("b_", term)
    if (!col %in% names(draws)) return(NULL)
    vals <- draws[[col]]
    tibble(
      term   = term,
      mean   = mean(vals),
      sd     = sd(vals),
      q0.025 = quantile(vals, 0.025),
      q0.975 = quantile(vals, 0.975),
      p_pos  = mean(vals > 0),
      sig_95 = as.integer(quantile(vals, 0.025) > 0 |
                            quantile(vals, 0.975) < 0)
    )
  }) |> bind_rows()
}

# ── Helper: pool posteriors across seeds ──────────────────────────────────────
pool_posteriors <- function(seed_rows) {
  seed_rows |>
    group_by(gene, response, term) |>
    summarise(
      n_seeds         = n(),
      mean_pooled     = mean(mean),
      sd_within       = mean(sd),
      sd_between      = sd(mean),
      sd_pooled       = sqrt(mean(sd^2) + var(mean)),
      p_pos_mean      = mean(p_pos),
      sig_95_majority = sum(sig_95),
      unstable        = as.integer(sd(mean) > mean(sd)),
      .groups = "drop"
    )
}

# ── Main loop ─────────────────────────────────────────────────────────────────
all_seed_rows <- list()
t_start_total <- proc.time()

for (gene in genes) {
  message("\n[", gene, "] brms multiseed analysis (", length(SEEDS),
          " seeds × n=", N_SUB, ")...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_out,   paste0(gene, "_brms_multiseed.csv"))

  if (!file.exists(rates_path)) { warning("Missing: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Missing: ", tree_path);  next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |>
    filter(accession %in% common, is.finite(log_rate)) |>
    arrange(match(accession, tree$tip.label))

  has_pop <- "scaled_pop_density" %in% names(df)
  if (has_pop) df <- df |> filter(!is.na(scaled_pop_density))

  tree <- keep.tip(tree, df$accession)

  predictors <- c("scaled_temp", "scaled_humidity",
                  if ("scaled_delta_temp"       %in% names(df)) "scaled_delta_temp",
                  if ("scaled_delta_humidity"   %in% names(df)) "scaled_delta_humidity",
                  if ("scaled_sampling_density" %in% names(df)) "scaled_sampling_density",
                  if (has_pop) "scaled_pop_density")

  message("  ", nrow(df), " taxa available")
  df <- add_phylo_clusters(df, tree)

  gene_seed_rows <- list()

  for (seed in SEEDS) {
    message("\n  [seed=", seed, "]")

    sub_df   <- stratified_subsample(df, N_SUB, seed = seed)
    sub_tree <- keep.tip(tree, sub_df$accession)
    sub_df   <- sub_df |> arrange(match(accession, sub_tree$tip.label))

    for (response in c("log_rate", "log_omega")) {
      if (!response %in% names(sub_df) ||
          all(is.na(sub_df[[response]]))) next

      sub_df_r   <- sub_df |> filter(is.finite(.data[[response]]))
      sub_tree_r <- keep.tip(sub_tree, sub_df_r$accession)
      A_r        <- vcv(sub_tree_r)
      A_r        <- A_r[sub_df_r$accession, sub_df_r$accession]
      A_r        <- A_r / max(diag(A_r))

      year_term <- if ("year" %in% names(sub_df_r)) "as.factor(year)" else NULL

      formula_brms <- as.formula(paste(
        response, "~",
        paste(c(predictors, year_term), collapse = " + "),
        "+ (1 | gr(accession, cov = A_r))"
      ))

      message("    Fitting brms: ", response,
              "  (n=", nrow(sub_df_r), ", ",
              CHAINS, " chains × ", ITER, " iter)")
      t0 <- proc.time()

      fit <- tryCatch(
        brm(
          formula  = formula_brms,
          data     = sub_df_r,
          data2    = list(A_r = A_r),
          family   = gaussian(),
          chains   = CHAINS,
          iter     = ITER,
          warmup   = ITER %/% 2,
          cores    = min(CHAINS, 16L),
          seed     = seed,
          silent   = 2,
          refresh  = 0,
          backend  = "cmdstanr",
          control  = list(adapt_delta = 0.95)
        ),
        error = function(e) {
          warning("brms failed for ", gene, "/", response,
                  "/seed=", seed, ": ", e$message)
          NULL
        }
      )

      if (is.null(fit)) next

      rhats  <- rhat(fit)
      bad_rh <- rhats[rhats > 1.01 & !is.na(rhats)]
      if (length(bad_rh) > 0)
        warning("R-hat > 1.01 for ", gene, "/", response, "/seed=", seed,
                ": ", paste(names(bad_rh), round(bad_rh, 3),
                            sep = "=", collapse = ", "))

      rows <- posterior_summary_brms(fit, predictors)
      rows <- bind_cols(
        tibble(gene = gene, response = response,
               seed = seed, n = nrow(sub_df_r)),
        rows,
        tibble(loo_elpd = tryCatch(
          loo(fit)$estimates["elpd_loo", "Estimate"],
          error = function(e) NA_real_
        ))
      )

      elapsed <- (proc.time() - t0)[["elapsed"]]
      rhat_max <- max(rhats[grepl("^b_", names(rhats))], na.rm = TRUE)
      message("    Done. Rhat max (fixed) = ", round(rhat_max, 3),
              "  [", round(elapsed / 60, 1), " min]")

      gene_seed_rows[[paste(seed, response, sep = "_")]] <- rows
    }
  }

  result_df <- bind_rows(gene_seed_rows)
  write_csv(result_df, out_path)
  message("\n  Written: ", out_path)
  all_seed_rows[[gene]] <- result_df
}

# ── Combine and pool ──────────────────────────────────────────────────────────
if (length(all_seed_rows) > 0) {
  combined_raw <- bind_rows(all_seed_rows)
  raw_path     <- file.path(paths_out, "all_brms_multiseed_raw.csv")
  write_csv(combined_raw, raw_path)
  message("\nWritten combined raw: ", raw_path)

  pooled      <- pool_posteriors(combined_raw)
  pooled_path <- file.path(paths_out, "all_brms_multiseed.csv")
  write_csv(pooled, pooled_path)
  message("Written pooled:       ", pooled_path)

  message("\n── Pooled posteriors (log_rate, predictor terms) ───────────────")
  pooled |>
    filter(term %in% PREDICTOR_TERMS, response == "log_rate") |>
    mutate(
      across(c(mean_pooled, sd_within, sd_between, sd_pooled),
             \(x) round(x, 4)),
      p_pos_mean = round(p_pos_mean, 3)
    ) |>
    arrange(gene, term) |>
    print(n = 60)

  n_unstable <- sum(pooled$unstable[pooled$term %in% PREDICTOR_TERMS],
                    na.rm = TRUE)
  message("\nUnstable predictors (sd_between > sd_within): ", n_unstable)
}

total_elapsed <- (proc.time() - t_start_total)[["elapsed"]]
message("\nTotal elapsed: ", round(total_elapsed / 60, 1), " min")
message("Done.")
