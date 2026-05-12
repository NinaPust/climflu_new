# scripts/03_rates_analysis/04_brms_phylo.R
#
# PURPOSE: Bayesian phylogenetic regression via brms (Stan MCMC).
#          Equivalent to 02_phylolm.R but returns full posteriors:
#          mean, SD, 95% credible interval, and exceedance probability
#          P(β > 0) for each predictor.
#
# APPROACH:
#   Trees of 9k–15k tips make dense VCV operations impractical, so tips are
#   stratified-subsampled by year (N_SUB = 2000). The phylogenetic covariance
#   matrix (vcv()) is passed to brms as a group-level correlation structure
#   via `gr(accession, cov = A)`, which implements a phylogenetic random
#   effect equivalent to Pagel's lambda. brms uses a non-centered
#   parameterization internally, which gives good mixing at this scale.
#
# RUNTIME: ~10–30 min per gene × response on 4 chains × 2000 iterations.
#   Set CHAINS and ITER below, or pass --chains and --iter on the command line.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_brms_results.csv   — posterior summary per coefficient
#   data/rates/all_brms_results.csv      — combined across all genes
#
# FLAGS:
#   --gene X       run a single gene segment
#   --test         use subsampled data from data/test/
#   --chains N     number of MCMC chains (default 4)
#   --iter N       total iterations per chain (default 2000)
#   --seed N       RNG seed for subsampling and Stan (default 42)

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
N_SUB <- Inf

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
SEED        <- as.integer(get_arg("--seed", 42))
seed_default <- 42L
seed_suffix  <- if (SEED != seed_default) paste0("_seed", SEED) else ""

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

# ── Helper: phylogenetically-stratified subsample ────────────────────────────
# Clusters tips by phylogenetic distance (cophenetic) using k-medoids, then
# stratifies the subsample by year × cluster so that both temporal coverage
# and tree topology are represented. Falls back to year-only stratification
# if the tree is not provided or n >= nrow(df).
#
# N_CLUST: number of phylogenetic clusters. More clusters = finer resolution
# but slower (cophenetic distance matrix is O(n^2)). 40 works well for
# influenza trees of 5k-10k tips.
stratified_subsample <- function(df, n, tree = NULL, seed = SEED,
                                 n_clust = 40) {
  set.seed(seed)
  if (nrow(df) <= n) return(df)

  # ── Phylogenetic clustering ───────────────────────────────────────────────
  if (!is.null(tree) && requireNamespace("cluster", quietly = TRUE)) {
    message("  Computing cophenetic distances for phylogenetic clustering...")
    coph   <- cophenetic(tree)
    # keep only tips in df (tree may be larger at this point)
    tips   <- intersect(rownames(coph), df$accession)
    coph   <- as.dist(coph[tips, tips])

    k      <- min(n_clust, length(tips) - 1L)
    pam_fit <- cluster::pam(coph, k = k, diss = TRUE)
    df$phylo_cluster <- pam_fit$clustering[df$accession]
    strat_vars <- c("phylo_cluster",
                    if ("year" %in% names(df)) "year")
  } else {
    strat_vars <- if ("year" %in% names(df)) "year" else NULL
  }

  # ── Proportional allocation across strata ────────────────────────────────
  if (is.null(strat_vars)) return(df[sample(nrow(df), n), ])

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

  # trim or top-up to exactly n
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

# ── Main loop ─────────────────────────────────────────────────────────────────
all_results <- list()
t_start_total <- proc.time()

for (gene in genes) {
  message("\n[", gene, "] brms phylogenetic regression...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_out, paste0(gene, "_brms_results", seed_suffix, ".csv"))

  if (!file.exists(rates_path)) { warning("Missing: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Missing: ", tree_path);  next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  # Align and clean (mirrors 02_phylolm.R)
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

  message("  ", nrow(df), " taxa available (using full dataset)")

  # sub_df   <- stratified_subsample(df, N_SUB, tree = tree)
  # sub_tree <- keep.tip(tree, sub_df$accession)
  # sub_df   <- sub_df |> arrange(match(accession, sub_tree$tip.label))
  sub_df   <- df
  sub_tree <- tree

  # Phylogenetic covariance matrix, scaled so max diagonal = 1
  A <- vcv(sub_tree)
  A <- A[sub_df$accession, sub_df$accession]
  A <- A / max(diag(A))

  gene_results <- list()

  for (response in c("log_rate", "log_omega")) {
    if (!response %in% names(sub_df) || all(is.na(sub_df[[response]]))) next
    sub_df_r  <- sub_df |> filter(is.finite(.data[[response]]))
    sub_tree_r <- keep.tip(sub_tree, sub_df_r$accession)
    A_r        <- vcv(sub_tree_r)
    A_r        <- A_r[sub_df_r$accession, sub_df_r$accession]
    A_r        <- A_r / max(diag(A_r))

    year_term <- if ("year" %in% names(sub_df)) "as.factor(year)" else NULL

    # IID residual variance is captured by sigma; a separate iid random effect
    # would be collinear with the phylogenetic term and cause non-identifiability.
    formula_brms <- as.formula(paste(
      response, "~",
      paste(c(predictors, year_term), collapse = " + "),
      "+ (1 | gr(accession, cov = A_r))"
    ))

    message("  Fitting brms: ", response, " ~ ", gene,
            "  (", CHAINS, " chains x ", ITER, " iter,",
            " n=", nrow(sub_df_r), ")")
    t_start <- proc.time()

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
        seed     = SEED,
        silent   = 2,
        refresh  = 0,
        backend  = "cmdstanr",
        control  = list(adapt_delta = 0.95)
      ),
      error = function(e) {
        warning("brms failed for ", gene, " / ", response, ": ", e$message)
        NULL
      }
    )

    if (is.null(fit)) next

    rhats  <- rhat(fit)
    bad_rh <- rhats[rhats > 1.01 & !is.na(rhats)]
    if (length(bad_rh) > 0)
      warning("R-hat > 1.01 for ", gene, "/", response, ": ",
              paste(names(bad_rh), round(bad_rh, 3),
                    sep = "=", collapse = ", "))

    rows <- posterior_summary_brms(fit, predictors)
    rows <- bind_cols(
      tibble(gene = gene, response = response, n = nrow(sub_df_r)),
      rows,
      tibble(loo_elpd = tryCatch(
        loo(fit)$estimates["elpd_loo", "Estimate"],
        error = function(e) NA_real_
      ))
    )

    gene_results[[response]] <- rows
    elapsed <- (proc.time() - t_start)[["elapsed"]]
    rhat_max <- max(rhats[grepl("^b_", names(rhats))], na.rm = TRUE)
    message("  Done. Rhat max (fixed) = ", round(rhat_max, 3),
            "  [", round(elapsed / 60, 1), " min]")
  }

  result_df <- bind_rows(gene_results)
  write_csv(result_df, out_path)
  message("  Written: ", out_path)

  result_df |>
    filter(term %in% PREDICTOR_TERMS) |>
    mutate(across(c(mean, sd, q0.025, q0.975), \(x) round(x, 4)),
           p_pos = round(p_pos, 3)) |>
    print()

  all_results[[gene]] <- result_df
}

# ── Combine ───────────────────────────────────────────────────────────────────
if (length(all_results) > 0) {
  combined <- bind_rows(all_results)
  combined_path <- file.path(paths_out,
                             paste0("all_brms_results", seed_suffix, ".csv"))
  write_csv(combined, combined_path)
  message("\nWritten combined: ", combined_path)

  message("\n── Predictor posteriors (log_rate, all genes) ──────────────────")
  combined |>
    filter(term %in% PREDICTOR_TERMS, response == "log_rate") |>
    mutate(
      ci    = paste0("[", round(q0.025, 3), ", ", round(q0.975, 3), "]"),
      p_pos = round(p_pos, 3),
      sig   = ifelse(sig_95 == 1, "*", "")
    ) |>
    select(gene, term, mean, ci, p_pos, sig) |>
    arrange(gene, term) |>
    print(n = 50)
}

total_elapsed <- (proc.time() - t_start_total)[["elapsed"]]
message("\nTotal elapsed: ", round(total_elapsed / 60, 1), " min")
message("Done.")
