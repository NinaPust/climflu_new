# scripts/03_rates_analysis/04b_inla_multiseed.R
#
# PURPOSE: Multi-seed variant of 04_inla_phylo.R.
#          Runs INLA over SEEDS × genes, then summarises:
#            - mean posterior (average across seeds)
#            - within-seed SD (average posterior SD)
#            - between-seed SD (SD of posterior means across seeds)
#          between-seed SD >> within-seed SD flags subsample-sensitive genes.
#
# APPROACH:
#   For each gene × seed, a fresh stratified subsample of N_SUB tips is drawn,
#   INLA is fitted, and the per-term posterior summary is stored. After all
#   seeds, summaries are pooled via moment-matching (mixture of Gaussians).
#   Only log_rate is modelled (log_omega skipped for speed).
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_inla_multiseed.csv   — per-seed posteriors
#   data/rates/all_inla_multiseed.csv      — pooled summary across seeds × genes
#
# FLAGS:
#   --gene X   run a single gene segment
#   --test     use subsampled data from data/test/

suppressPackageStartupMessages({
  if (!requireNamespace("INLA", quietly = TRUE))
    stop(
      "INLA is not installed. Run this in R first:\n\n",
      "  install.packages(\"INLA\",\n",
      "    repos = c(getOption(\"repos\"),\n",
      "              INLA = \"https://inla.r-inla-download.org/R/stable\"),\n",
      "    dep = TRUE)\n"
    )
  library(INLA)
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
  library(Matrix)
})

source(here("config.R"))
source(here("R", "climate_utils.R"))

# ── Settings ──────────────────────────────────────────────────────────────────
N_SUB <- 1000
SEEDS <- c(42, 123, 7, 999, 2024, 314, 888)

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args      <- commandArgs(trailingOnly = TRUE)
genes     <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}
use_test <- "--test" %in% args

if (use_test) message("*** Running in TEST mode ***")

paths_timed <- if (use_test) PATHS$test_trees_timed else PATHS$trees_timed
paths_rates <- if (use_test) PATHS$test_rates       else PATHS$rates

PREDICTOR_TERMS <- c("scaled_temp", "scaled_humidity",
                     "scaled_delta_temp", "scaled_delta_humidity",
                     "scaled_sampling_density", "scaled_pop_density")

# ── Helper: extract INLA posterior summary for one term ───────────────────────
posterior_summary <- function(fit, term) {
  marg <- fit$marginals.fixed[[term]]
  if (is.null(marg)) return(NULL)

  q <- inla.qmarginal(c(0.025, 0.975), marg)
  tibble(
    mean   = inla.emarginal(function(x) x, marg),
    sd     = sqrt(inla.emarginal(function(x) x^2, marg) -
                    inla.emarginal(function(x) x, marg)^2),
    q0.025 = q[1],
    q0.975 = q[2],
    p_pos  = 1 - inla.pmarginal(0, marg),
    sig_95 = as.integer(q[1] > 0 | q[2] < 0)
  )
}

# ── Helper: stratified subsample (by year) ────────────────────────────────────
stratified_subsample <- function(df, n, seed) {
  set.seed(seed)
  if (nrow(df) <= n) return(df)
  total_n     <- nrow(df)
  year_counts <- df |>
    count(.data$year, name = "grp_n") |>
    mutate(.k = pmax(1L, as.integer(round(.data$grp_n / total_n * n))))
  sampled <- lapply(split(df, df$year), function(sub) {
    k <- year_counts$.k[year_counts$year == sub$year[1]]
    sub[sample(nrow(sub), min(k, nrow(sub))), ]
  })
  result <- do.call(rbind, sampled)
  result[sample(nrow(result), min(n, nrow(result))), ]
}

# ── Pool posteriors across seeds (moment-matching) ───────────────────────────
# For a mixture of K Gaussians with means μ_k and SDs σ_k:
#   pooled mean  = mean(μ_k)
#   pooled SD    = sqrt(mean(σ_k²) + var(μ_k))   [law of total variance]
#   between-seed SD = sd(μ_k)   — key instability diagnostic
pool_posteriors <- function(seed_rows) {
  seed_rows |>
    group_by(.data$gene, .data$response, .data$term) |>
    summarise(
      n_seeds        = n(),
      n              = mean(.data$n),
      # pooled posterior mean
      mean           = mean(.data$mean),
      # pooled posterior SD (law of total variance)
      sd_pooled      = sqrt(mean(.data$sd^2) + var(.data$mean)),
      # within-seed SD (average per-seed uncertainty)
      sd_within      = mean(.data$sd),
      # between-seed SD (subsample sensitivity)
      sd_between     = sd(.data$mean),
      # pooled exceedance prob (average across seeds)
      p_pos          = mean(.data$p_pos),
      # significant in majority of seeds
      sig_95_votes   = sum(.data$sig_95),
      sig_95_majority = as.integer(sum(.data$sig_95) > n() / 2),
      # mean DIC/WAIC for reference
      dic            = mean(.data$dic),
      waic           = mean(.data$waic),
      .groups = "drop"
    ) |>
    mutate(
      # 95% CrI from pooled mean ± 1.96 × pooled SD
      q0.025 = .data$mean - 1.96 * .data$sd_pooled,
      q0.975 = .data$mean + 1.96 * .data$sd_pooled,
      # flag if between-seed variance dominates within-seed variance
      unstable = as.integer(.data$sd_between > .data$sd_within)
    )
}

# ── Main loop ─────────────────────────────────────────────────────────────────
all_seed_rows <- list()

for (gene in genes) {
  message("\n[", gene, "] Multi-seed INLA (", length(SEEDS), " seeds × n=", N_SUB, ")...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_rates, paste0(gene, "_inla_multiseed.csv"))

  if (!file.exists(rates_path)) { warning("Missing: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Missing: ", tree_path);  next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  # Clean (mirrors 04_inla_phylo.R)
  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |>
    filter(.data$accession %in% common, is.finite(.data$log_rate)) |>
    arrange(match(.data$accession, tree$tip.label))

  has_pop <- "scaled_pop_density" %in% names(df)
  if (has_pop) df <- df |> filter(!is.na(.data$scaled_pop_density))

  tree <- keep.tip(tree, df$accession)

  predictors <- c("scaled_temp", "scaled_humidity",
                  if ("scaled_delta_temp"       %in% names(df)) "scaled_delta_temp",
                  if ("scaled_delta_humidity"   %in% names(df)) "scaled_delta_humidity",
                  if ("scaled_sampling_density" %in% names(df)) "scaled_sampling_density",
                  if (has_pop) "scaled_pop_density")

  message("  ", nrow(df), " taxa available after cleaning")

  gene_seed_rows <- list()

  for (seed in SEEDS) {
    message("  seed=", seed, " ...")

    sub_df   <- stratified_subsample(df, N_SUB, seed)
    sub_tree <- keep.tip(tree, sub_df$accession)
    sub_df   <- sub_df |> arrange(match(.data$accession, sub_tree$tip.label))

    V <- vcv(sub_tree)
    V <- V[sub_df$accession, sub_df$accession]
    V <- V / max(diag(V))
    P <- Matrix(solve(V), sparse = FALSE)

    sub_df$tip_id_phylo <- seq_len(nrow(sub_df))
    sub_df$tip_id_iid   <- seq_len(nrow(sub_df))
    rownames(sub_df)    <- sub_df$accession

    year_term <- if ("year" %in% names(sub_df)) "as.factor(year)" else NULL

    formula_inla <- as.formula(paste(
      "log_rate ~",
      paste(c(predictors, year_term), collapse = " + "),
      "+ f(tip_id_phylo, model = 'generic0', Cmatrix = P,",
      "    hyper = list(prec = list(prior = 'pc.prec',",
      "                             param = c(1, 0.01))))",
      "+ f(tip_id_iid, model = 'iid',",
      "    hyper = list(prec = list(prior = 'pc.prec',",
      "                             param = c(1, 0.01))))"
    ))

    fit <- tryCatch(
      inla(
        formula_inla,
        data    = sub_df,
        family  = "gaussian",
        control.compute = list(dic = TRUE, waic = TRUE, cpo = FALSE),
        control.predictor = list(compute = FALSE),
        verbose = FALSE
      ),
      error = function(e) {
        warning("INLA failed for ", gene, " seed=", seed, ": ", e$message)
        NULL
      }
    )

    if (is.null(fit)) next

    all_terms <- c("(Intercept)", predictors)
    rows <- lapply(all_terms, function(term) {
      s <- posterior_summary(fit, term)
      if (is.null(s)) return(NULL)
      bind_cols(
        tibble(gene = gene, response = "log_rate", term = term,
               seed = seed, n = nrow(sub_df)),
        s,
        tibble(dic = fit$dic$dic, waic = fit$waic$waic)
      )
    })

    gene_seed_rows[[as.character(seed)]] <- bind_rows(rows)
    message("    DIC=", round(fit$dic$dic, 1),
            "  WAIC=", round(fit$waic$waic, 1))
  }

  if (length(gene_seed_rows) == 0) next

  gene_df <- bind_rows(gene_seed_rows)
  write_csv(gene_df, out_path)
  message("  Written per-seed: ", out_path)

  all_seed_rows[[gene]] <- gene_df
}

# ── Pool and summarise ────────────────────────────────────────────────────────
if (length(all_seed_rows) == 0) {
  message("No results to pool.")
  quit(save = "no")
}

all_raw      <- bind_rows(all_seed_rows)
combined_raw <- file.path(paths_rates, "all_inla_multiseed_raw.csv")
write_csv(all_raw, combined_raw)
message("\nWritten raw (per-seed): ", combined_raw)

pooled      <- pool_posteriors(all_raw)
combined_pooled <- file.path(paths_rates, "all_inla_multiseed.csv")
write_csv(pooled, combined_pooled)
message("Written pooled: ", combined_pooled)

message("\n── Pooled predictor posteriors (log_rate) ───────────────────────")
pooled |>
  filter(.data$term %in% PREDICTOR_TERMS, .data$response == "log_rate") |>
  mutate(
    ci         = paste0("[", round(.data$q0.025, 3), ", ", round(.data$q0.975, 3), "]"),
    p_pos      = round(.data$p_pos, 3),
    sd_btwn    = round(.data$sd_between, 4),
    sig        = ifelse(.data$sig_95_majority == 1, "*", ""),
    unstable   = ifelse(.data$unstable == 1, "!", "")
  ) |>
  select(.data$gene, .data$term, .data$mean, .data$ci, .data$p_pos,
         .data$sig, .data$sd_btwn, .data$unstable, .data$n_seeds) |>
  arrange(.data$gene, .data$term) |>
  print(n = 60)

message("\nNote: '!' in unstable column = between-seed SD > within-seed SD",
        " (subsample-sensitive estimate)")
message("\nDone.")
