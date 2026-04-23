# scripts/03_rates_analysis/04_inla_phylo.R
#
# PURPOSE: Bayesian phylogenetic regression via INLA.
#          Equivalent to 02_phylolm.R but returns full posterior distributions
#          (mean, SD, 95% credible interval, exceedance probability P(β > 0))
#          rather than point estimates and p-values.
#
# APPROACH:
#   For each gene, N_SUB tips are randomly sampled (stratified by year so
#   temporal coverage is preserved). The phylogenetic VCV matrix is computed
#   on this subsample and inverted to a precision matrix, which is passed to
#   INLA as a generic0 random effect. A second iid random effect absorbs
#   non-phylogenetic residual variance. Together these two random effects
#   are equivalent to Pagel's lambda decomposition.
#
#   Tree sizes (9k–15k tips) make full-data dense matrix operations impractical
#   (~800 MB VCV). N_SUB = 1000 gives manageable 1 MB matrices while
#   retaining good statistical power.
#
# INSTALL INLA (run once in R before using this script):
#   install.packages("INLA",
#     repos = c(getOption("repos"),
#               INLA = "https://inla.r-inla-download.org/R/stable"),
#     dep = TRUE)
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_inla_results.csv  — posterior summary per coefficient
#   data/rates/all_inla_results.csv     — combined across all genes
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

# ── Settings ───────────────────────────────────────────────────────────────────
# Number of tips to subsample per gene. Increase for more power (slower);
# decrease if memory is tight. 1000 → ~8 MB VCV matrix.
N_SUB <- 1000

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
                     "scaled_pop_density", "scaled_season")

# ── Helper: extract INLA posterior summary for one term ───────────────────────
posterior_summary <- function(fit, term) {
  marg <- fit$marginals.fixed[[term]]
  if (is.null(marg)) return(NULL)

  q <- inla.qmarginal(c(0.025, 0.975), marg)
  tibble(
    mean      = inla.emarginal(function(x) x, marg),
    sd        = sqrt(inla.emarginal(function(x) x^2, marg) -
                       inla.emarginal(function(x) x, marg)^2),
    q0.025    = q[1],
    q0.975    = q[2],
    # P(β > 0): directional evidence — Bayesian analogue of one-sided p-value
    p_pos     = 1 - inla.pmarginal(0, marg),
    # 1 if 95% CrI excludes 0, 0 otherwise
    sig_95    = as.integer(q[1] > 0 | q[2] < 0)
  )
}

# ── Helper: stratified subsample (by year) ────────────────────────────────────
stratified_subsample <- function(df, n, seed = 42) {
  set.seed(seed)
  if (nrow(df) <= n) return(df)
  # Compute per-year sample sizes up front, then sample each group by integer
  total_n <- nrow(df)
  year_counts <- df |>
    count(year, name = "grp_n") |>
    mutate(.k = pmax(1L, as.integer(round(grp_n / total_n * n))))

  sampled <- lapply(split(df, df$year), function(sub) {
    k <- year_counts$.k[year_counts$year == sub$year[1]]
    sub[sample(nrow(sub), min(k, nrow(sub))), ]
  })
  result <- do.call(rbind, sampled)
  # Top up / trim to exactly n
  result[sample(nrow(result), min(n, nrow(result))), ]
}

# ── Main loop ─────────────────────────────────────────────────────────────────
all_results <- list()

for (gene in genes) {
  message("\n[", gene, "] INLA phylogenetic regression...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_rates, paste0(gene, "_inla_results.csv"))

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
                  if (has_pop) "scaled_pop_density",
                  if ("scaled_season" %in% names(df)) "scaled_season")

  message("  ", nrow(df), " taxa available; subsampling to ",
          min(N_SUB, nrow(df)))

  # ── Subsample ───────────────────────────────────────────────────────────────
  sub_df   <- stratified_subsample(df, N_SUB)
  sub_tree <- keep.tip(tree, sub_df$accession)
  sub_df   <- sub_df |>
    arrange(match(accession, sub_tree$tip.label))

  # ── Phylogenetic precision matrix ────────────────────────────────────────────
  # vcv() returns the Brownian-motion covariance matrix (shared branch lengths)
  V <- vcv(sub_tree)
  V <- V[sub_df$accession, sub_df$accession]   # reorder to match data
  V <- V / max(diag(V))                         # scale so max variance = 1

  # Invert to precision matrix; wrap in Matrix class for INLA
  P <- Matrix(solve(V), sparse = FALSE)

  # Integer tip index (1-based, matching rows of sub_df / V)
  sub_df$tip_id_phylo <- seq_len(nrow(sub_df))
  sub_df$tip_id_iid   <- seq_len(nrow(sub_df))

  rownames(sub_df) <- sub_df$accession

  gene_results <- list()

  for (response in c("log_rate", "log_omega")) {
    if (!response %in% names(sub_df) ||
        all(is.na(sub_df[[response]]))) next

    formula_inla <- as.formula(paste(
      response, "~",
      paste(predictors, collapse = " + "),
      # Phylogenetic random effect — covariance proportional to VCV
      "+ f(tip_id_phylo, model = 'generic0', Cmatrix = P,",
      "    hyper = list(prec = list(prior = 'pc.prec',",
      "                             param = c(1, 0.01))))",
      # IID residual — absorbs non-phylogenetic variance
      "+ f(tip_id_iid, model = 'iid',",
      "    hyper = list(prec = list(prior = 'pc.prec',",
      "                             param = c(1, 0.01))))"
    ))

    message("  Fitting INLA: ", response, " ~ ", gene)

    fit <- tryCatch(
      inla(
        formula_inla,
        data    = sub_df,
        family  = "gaussian",
        control.compute = list(
          dic  = TRUE,
          waic = TRUE,
          cpo  = FALSE   # expensive; enable if needed
        ),
        control.predictor = list(compute = FALSE),
        verbose = FALSE
      ),
      error = function(e) {
        warning("INLA failed for ", gene, " / ", response, ": ", e$message)
        NULL
      }
    )

    if (is.null(fit)) next

    # Extract posterior summaries for all terms
    all_terms <- c("(Intercept)", predictors)
    rows <- lapply(all_terms, function(term) {
      s <- posterior_summary(fit, term)
      if (is.null(s)) return(NULL)
      bind_cols(
        tibble(gene = gene, response = response, term = term, n = nrow(sub_df)),
        s,
        tibble(dic = fit$dic$dic, waic = fit$waic$waic)
      )
    })

    gene_results[[response]] <- bind_rows(rows)
    message("  Done. DIC = ", round(fit$dic$dic, 1),
            "  WAIC = ", round(fit$waic$waic, 1))
  }

  result_df <- bind_rows(gene_results)
  write_csv(result_df, out_path)
  message("  Written: ", out_path)

  # Print predictor terms only
  result_df |>
    filter(term %in% PREDICTOR_TERMS) |>
    mutate(across(c(mean, sd, q0.025, q0.975), \(x) round(x, 4)),
           p_pos = round(p_pos, 3)) |>
    print()

  all_results[[gene]] <- result_df
}

# ── Combine and summarise ─────────────────────────────────────────────────────
if (length(all_results) > 0) {
  combined <- bind_rows(all_results)
  combined_path <- file.path(paths_rates, "all_inla_results.csv")
  write_csv(combined, combined_path)
  message("\nWritten combined: ", combined_path)

  message("\n── Predictor posteriors (log_rate, all genes) ──────────────────")
  combined |>
    filter(term %in% PREDICTOR_TERMS, response == "log_rate") |>
    mutate(
      ci       = paste0("[", round(q0.025, 3), ", ", round(q0.975, 3), "]"),
      p_pos    = round(p_pos, 3),
      sig      = ifelse(sig_95 == 1, "*", "")
    ) |>
    select(gene, term, mean, ci, p_pos, sig) |>
    arrange(gene, term) |>
    print(n = 50)
}

message("\nDone.")
