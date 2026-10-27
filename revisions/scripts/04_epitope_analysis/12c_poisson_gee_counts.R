# revisions/scripts/04_epitope_analysis/12c_poisson_gee_counts.R
#
# SUPPORTING ANALYSIS — count model that does NOT split the hurdle.
#
# PURPOSE: Poisson regression of non-synonymous substitution counts against
#          climate predictors, with log-exposure as a covariate. Third leg of
#          the revision, with a specific job:
#
#            12a  binary "any change"        keeps zeros, loses magnitude
#            12b  binary "was it NS"         keeps type, drops zero-change tips
#            12c  count of NS substitutions  keeps BOTH
#
#          Because it keeps both, it checks that the two-part structure is not
#          manufacturing the result. If 12a and 12b agree in direction with 12c,
#          the hurdle is a presentational split of one coherent signal. If they
#          disagree, the split is doing something that needs explaining.
#
# WHY TWO MODELS HERE:
#          phylolm's poisson_GEE does not scale to this dataset. Measured on
#          H1 whole gene (identical formula, 4 BLAS threads):
#
#            n = 1500   converges, exposure = 0.199
#            n = 3000   converges, exposure = 0.250
#            n = 5000   "system is computationally singular"
#            n = 9311   convergence code 1, exposure = -254, SE = NaN
#
#          It breaks between 3k and 5k tips. Centring the exposure and dropping
#          the year factor do not help, and the design matrix is well
#          conditioned (kappa = 235), so this is not collinearity in the
#          predictors. With mean counts of 0.46 (whole gene) down to 0.011
#          (ridge), the fitted means span so many orders of magnitude that the
#          GEE's iterative working matrix goes numerically singular. phyloglm
#          exposes no tolerance argument for that path (btol is logistic-only).
#
#          So this script fits two things:
#
#            GLM_full_n   Poisson GLM, ALL tips, NO phylogenetic correction.
#                         Fast and stable (dispersion ~0.97, so the Poisson
#                         assumption itself is sound). This is the primary
#                         result and the one plotted.
#
#            GEE_subN     poisson_GEE on --gee-reps random subsamples of
#                         --gee-n tips, WITH phylogenetic correction. Reported
#                         as the spread of coefficients across replicates. Its
#                         job is narrow: show that omitting the phylogenetic
#                         correction in GLM_full_n does not flip signs.
#
#          Note the two are not strictly comparable: subsampling reassigns
#          terminal branch lengths in the pruned tree and changes the
#          correlation structure, so treat GEE_subN as a direction/sign check,
#          not as a second estimate of the same quantity.
#
# LIMITATION TO STATE IN THE PAPER: the primary 12c fit does not correct for
#          phylogeny. That is acceptable for its role — 12a and 12b both carry
#          phylogenetic correction, and 12c exists to test the hurdle split
#          rather than to provide an independent effect estimate.
#
# EXPOSURE: log(N_sites * branch_time) — non-synonymous sites only, because the
#          response counts only non-synonymous substitutions. (Contrast 12a,
#          whose response includes synonymous changes and so uses total sites.)
#          Under a clean Poisson process its coefficient should be near 1;
#          observed values are well below that (~0.46 whole gene), consistent
#          with 12a and with the earlier 12d check, so it is printed per fit as
#          a diagnostic about the process rather than treated as a nuisance.
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/phylolm/{gene}_poisson_nscount_results.csv
#   data/rates/phylolm/all_poisson_nscount_results.csv
#   data/rates/phylolm/all_poisson_gee_subsample_check.csv
#   figures/epitope/poisson_nscount_forest.png
#   figures/epitope/poisson_nscount_forest_fdr.png
#
# FLAGS:
#   --gene X       run a single gene (default: H1 H3)
#   --region X     run a single region (default: full HA2 epitope ridge)
#   --no-year      drop as.factor(year)
#   --with-gee     ALSO run the phylogenetic GEE subsample check. OFF BY
#                  DEFAULT: at n = 3000 only about 1 replicate in 4 converges,
#                  so the check is too thin to be worth its runtime as it
#                  stands. The code is kept because a smaller --gee-n may
#                  converge reliably enough to be useful; that has not been
#                  established, so nothing in the paper should rest on it yet.
#   --gee-n N      tips per GEE replicate (default 3000; 5000 is known to fail)
#   --gee-reps N   number of GEE replicates (default 10)
#   --threads N    BLAS/OMP cap (default 4; shared machine)
#   --seed N       default 42
#
# RUNTIME: the GLMs are near-instant. The GEE replicates are the only real cost
#   (~2 s each at n = 3000), so the whole script runs in minutes, not hours.

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "phylolm_utils.R"))
source(here("revisions", "R", "hurdle_utils.R"))
source(here("revisions", "R", "forest_utils.R"))

args     <- commandArgs(trailingOnly = TRUE)
genes    <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
regions  <- if ("--region" %in% args) args[which(args == "--region") + 1] else
            c("full", "HA2", "epitope", "ridge")
use_year <- !("--no-year"  %in% args)
skip_gee <- "--skip-gee"   %in% args
gee_n    <- if ("--gee-n"    %in% args) as.integer(args[which(args == "--gee-n")    + 1]) else 3000L
gee_reps <- if ("--gee-reps" %in% args) as.integer(args[which(args == "--gee-reps") + 1]) else 10L
n_thread <- if ("--threads"  %in% args) as.integer(args[which(args == "--threads")  + 1]) else 4L
seed     <- if ("--seed"     %in% args) as.integer(args[which(args == "--seed")     + 1]) else 42L

cap_blas_threads(n_thread)
set.seed(seed)

if (gee_n >= 5000)
  warning("--gee-n ", gee_n, " is at or above the size where poisson_GEE was ",
          "observed to fail (5000). Expect replicates to error out.",
          immediate. = TRUE)

predictor_terms <- c(
  "scaled_temp", "scaled_humidity",
  "scaled_delta_temp", "scaled_delta_humidity",
  "scaled_sampling_density", "scaled_pop_density"
)

out_dir <- file.path(PATHS$rates_phylolm)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


#' Tidy a Poisson GLM into the shared results schema.
#'
#' Schema matches run_phyloglm() so GLM, phyloglm and GEE results can be
#' stacked and plotted with the same code. alpha/logLik-style phylogenetic
#' fields are NA by construction: this fit has no phylogenetic component.
tidy_glm <- function(fit, gene, region, n_obs) {
  cc <- as.data.frame(summary(fit)$coefficients)
  names(cc)[1:4] <- c("estimate", "std_error", "z_value", "p.value")
  cc$term <- rownames(cc); rownames(cc) <- NULL
  cc$gene      <- gene
  cc$subset    <- region
  cc$response  <- "ns_count"
  cc$model     <- "GLM_full_n"
  cc$alpha     <- NA_real_
  cc$logLik    <- as.numeric(logLik(fit))
  cc$AIC       <- AIC(fit)
  cc$AIC_null  <- NA_real_
  cc$n         <- n_obs
  cc$n_events  <- sum(fit$y > 0)
  cc$dispersion <- sum(residuals(fit, "pearson")^2) / fit$df.residual
  cc$r.squared <- NA_real_
  cc[, c("gene", "subset", "response", "model", "term",
         "estimate", "std_error", "p.value", "alpha", "logLik",
         "AIC", "AIC_null", "n", "n_events", "dispersion", "r.squared")]
}

# ── Fit per gene x region ─────────────────────────────────────────────────────

gee_check <- list()

for (gene in genes) {
  message("\n[", gene, "] Poisson count models — NS substitutions")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
  if (!file.exists(rates_path)) { warning("Rates not found: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Tree not found: ",  tree_path);  next }

  df          <- read_csv(rates_path, show_col_types = FALSE)
  full_tree   <- read.tree(tree_path)
  branch_time <- terminal_branch_times(full_tree)

  gene_results <- list()

  for (r in regions) {
    ns_col <- paste0(r, "_ns_count"); N_col <- paste0(r, "_N_sites")
    if (!all(c(ns_col, N_col) %in% names(df))) {
      message("  [", r, "] missing columns — re-run 09_substitutions.R; skipping")
      next
    }

    sub <- df
    sub$ns_count     <- sub[[ns_col]]
    sub$log_exposure <- log(sub[[N_col]] * branch_time[sub$accession])

    keep <- intersect(c("accession", "ns_count", "log_exposure",
                        predictor_terms, if (use_year) "year"), names(sub))
    sub  <- sub[, keep, drop = FALSE]
    sub  <- sub[stats::complete.cases(sub) &
                is.finite(sub$log_exposure), , drop = FALSE]
    if (nrow(sub) < 50) { message("  [", r, "] <50 usable tips, skipping"); next }

    n_nonzero <- sum(sub$ns_count > 0)
    message("  [", r, "] n = ", nrow(sub), "  tips with NS > 0 = ", n_nonzero,
            " (", round(100 * n_nonzero / nrow(sub), 1), "%)",
            "  mean count = ", round(mean(sub$ns_count), 3))
    if (n_nonzero < 20) {
      message("    fewer than 20 tips with a non-synonymous change, skipping")
      next
    }

    preds     <- intersect(predictor_terms, names(sub))
    year_term <- if (use_year && "year" %in% names(sub)) "as.factor(year)" else NULL
    fml <- as.formula(paste("ns_count ~",
                            paste(c(preds, "log_exposure", year_term),
                                  collapse = " + ")))

    # ── Primary: Poisson GLM at full n ────────────────────────────────────────
    fit <- tryCatch(glm(fml, data = sub, family = poisson),
                    error = function(e) {
                      message("    GLM FAILED: ", conditionMessage(e)); NULL })
    if (is.null(fit)) next

    res <- tidy_glm(fit, gene, r, nrow(sub))
    disp <- res$dispersion[1]
    exp_est <- res$estimate[res$term == "log_exposure"]
    message("    GLM ok  |  exposure = ", signif(exp_est, 3),
            " (Poisson expectation ~1)  |  dispersion = ", signif(disp, 3))
    if (disp > 1.5)
      message("    NOTE: dispersion > 1.5 — counts are overdispersed for ",
              "Poisson; consider quasipoisson or negative binomial here")
    gene_results[[r]] <- res

    # ── Secondary: phylogenetic GEE on subsamples ─────────────────────────────
    # Narrow purpose: does adding phylogenetic correction flip any signs?
    if (skip_gee) next
    if (nrow(sub) <= gee_n) {
      message("    GEE check: n <= ", gee_n, ", fitting once on all tips")
      reps <- 1L
    } else {
      reps <- gee_reps
    }

    ok_reps <- 0L
    for (i in seq_len(reps)) {
      s0 <- if (nrow(sub) <= gee_n) sub else sub[sample(nrow(sub), gee_n), ]
      p  <- tryCatch(prep_tree_data(s0, full_tree, quiet = TRUE),
                     error = function(e) NULL)
      if (is.null(p)) next
      rr <- tryCatch(
        run_phyloglm_gee(fml, p$data, p$tree, subset = r,
                         response = "ns_count", gene = gene),
        error = function(e) NULL
      )
      if (is.null(rr)) next
      rr <- rr[rr$term %in% c(predictor_terms, "log_exposure"), ]
      rr$rep   <- i
      rr$gee_n <- nrow(p$data)
      gee_check[[paste(gene, r, i, sep = "_")]] <- rr
      ok_reps <- ok_reps + 1L
    }
    message("    GEE check: ", ok_reps, "/", reps, " replicates converged",
            if (ok_reps == 0) "  (region too sparse for a count model)" else "")
  }

  if (length(gene_results) == 0) { message("  no successful fits"); next }
  gene_df  <- bind_rows(gene_results)
  gene_out <- file.path(out_dir, paste0(gene, "_poisson_nscount_results.csv"))
  write_csv(gene_df, gene_out)
  message("  Written: ", basename(gene_out))
}

# ── Combine, FDR-correct, report ──────────────────────────────────────────────

all_res <- lapply(genes, function(g) {
  p <- file.path(out_dir, paste0(g, "_poisson_nscount_results.csv"))
  if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL
}) |> bind_rows()

if (nrow(all_res) == 0) stop("No results produced.")

all_res <- add_fdr(all_res, predictor_terms)
combined <- file.path(out_dir, "all_poisson_nscount_results.csv")
write_csv(all_res, combined)
message("\nWritten: ", basename(combined))

message("\nExposure coefficients (Poisson expectation ~1) and dispersion:")
all_res |>
  filter(term == "log_exposure") |>
  mutate(estimate = round(estimate, 3), dispersion = round(dispersion, 3)) |>
  select(gene, subset, estimate, std_error, dispersion, n) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nSignificant after BH-FDR (p.fdr < 0.05):")
sig <- all_res |>
  filter(term %in% predictor_terms, !is.na(p.fdr), p.fdr < 0.05) |>
  mutate(rate_multiplier = round(exp(estimate), 3),
         estimate = round(estimate, 4),
         across(c(p.value, p.fdr), \(x) signif(x, 3))) |>
  select(gene, subset, term, estimate, rate_multiplier, p.value, p.fdr, n) |>
  arrange(p.fdr)
if (nrow(sig) == 0) message("  (none)") else print(as.data.frame(sig))

# ── GEE subsample check: does phylogenetic correction flip any signs? ─────────

if (length(gee_check) > 0) {
  gee_df <- bind_rows(gee_check)
  gee_out <- file.path(out_dir, "all_poisson_gee_subsample_check.csv")
  write_csv(gee_df, gee_out)
  message("\nWritten: ", basename(gee_out))

  glm_signs <- all_res |>
    filter(term %in% c(predictor_terms, "log_exposure")) |>
    select(gene, subset, term, glm_estimate = estimate)

  message("\nGEE subsample check (", gee_reps, " reps of n = ", gee_n,
          ", phylogenetically corrected) vs GLM at full n:")
  gee_df |>
    group_by(gene, subset, term) |>
    summarise(reps = n(),
              gee_median = median(estimate),
              gee_lo = min(estimate), gee_hi = max(estimate),
              .groups = "drop") |>
    left_join(glm_signs, by = c("gene", "subset", "term")) |>
    mutate(sign_agrees = sign(gee_median) == sign(glm_estimate),
           across(c(gee_median, gee_lo, gee_hi, glm_estimate),
                  \(x) round(x, 4))) |>
    as.data.frame() |>
    print(row.names = FALSE)
} else if (!skip_gee) {
  message("\nGEE subsample check: no replicates converged in any region.")
}

# ── Forest plots ──────────────────────────────────────────────────────────────

save_forest_pair(
  all_res,
  out_dir = here("revisions", "figures", "epitope"),
  stem    = "poisson_nscount_forest",
  title   = "Supporting: Poisson count model of non-synonymous substitutions",
  xlab    = "Log-rate coefficient"
)

message("\nDone.")
