# scripts/03_rates_analysis/07_phylolm_jackknife.R
#
# PURPOSE: Block jackknife subsampling for phylolm temperature analysis.
#          Runs 100 random samples at four draw sizes (1000 fixed, 15%, 30%,
#          60% of dataset), runs phylolm on each, and saves per-gene results for
#          downstream robustness plots.
#
#          REVISED: fits the FULL model (six climate predictors + year factor),
#          matching 02_phylolm.R, instead of the submitted temperature-only
#          formula. Results are long format, one row per replicate per term.
#
# INPUT  (from the ORIGINAL data tree — the revision did not change these):
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT goes to data/rates/, leaving the submitted jackknife intact.
#
# OUTPUT (per draw size):
#   data/rates/jackknife/{gene}/jackknife_results.csv          # n=1000 (legacy path)
#   data/rates/jackknife/{gene}/jackknife_sequence_membership.csv
#   data/rates/jackknife/{gene}/pct15/jackknife_results.csv    # 15% of dataset
#   data/rates/jackknife/{gene}/pct30/jackknife_results.csv    # 30% of dataset
#   data/rates/jackknife/{gene}/pct60/jackknife_results.csv    # 60% of dataset
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
  library(tidyr)
})

source(here("config.R"))
source(here("R", "climate_utils.R"))

N_SAMPLES   <- 1000L
MASTER_SEED <- 42L

# REVISED: the full model, matching 03_rates_analysis/02_phylolm.R. The
# submitted version jackknifed `log_rate ~ scaled_temp + as.factor(year)`, a
# reduced formula, so its temperature coefficient did not mean the same thing as
# the one reported in the figures — with the other climate terms absent, the
# temperature coefficient absorbed whatever they would have explained. Every
# coefficient is now stored (long format, one row per replicate per term) so any
# predictor can be plotted later without refitting.
PRED <- c("scaled_temp", "scaled_humidity",
          "scaled_delta_temp", "scaled_delta_humidity",
          "scaled_sampling_density", "scaled_pop_density")

# Draw sizes: fixed count (n1000) or fraction of dataset (pct15/30/60).
# n1000 uses the legacy top-level path; percentage sizes use a subdirectory.
DRAW_SIZES <- list(
  list(label = "n1000", fixed = 1000L, frac = NA),
  list(label = "pct15", fixed = NA,    frac = 0.15),
  list(label = "pct30", fixed = NA,    frac = 0.30),
  list(label = "pct60", fixed = NA,    frac = 0.60)
)

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  GENES
}

for (gene in genes) {
  message("\n[", gene, "] Loading data...")

  rates_path <- here("data", "rates", paste0(gene, "_rates.csv"))
  tree_path  <- here("data", "trees", "timed", paste0(gene, "_timetree.nwk"))

  if (!file.exists(rates_path)) {
    warning("Rates file not found: ", rates_path)
    next
  }
  if (!file.exists(tree_path)) {
    warning("Tree not found: ", tree_path)
    next
  }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  df   <- scale_env_vars(df)
  tree <- read.tree(tree_path)

  # ── Align tips and data ────────────────────────────────────────────────────
  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  # Drop zero-length terminal branches (phylolm singularity)
  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    message("  Dropping ", length(zero_tips), " tips with zero-length branches")
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |>
    filter(accession %in% common) |>
    filter(!is.na(log_rate), if_all(all_of(intersect(PRED, names(df))), \(x) !is.na(x))) |>
    arrange(match(accession, tree$tip.label))
  tree <- keep.tip(tree, df$accession)
  rownames(df) <- df$accession

  n_avail <- nrow(df)
  preds      <- intersect(PRED, names(df))
  year_term  <- if ("year" %in% names(df)) "as.factor(year)" else NULL
  jk_formula <- as.formula(paste("log_rate ~",
                                 paste(c(preds, year_term), collapse = " + ")))
  message("  ", n_avail, " taxa available for subsampling")
  message("  formula: ", deparse(jk_formula))

  if (n_avail < 10) {
    warning("  Too few taxa for jackknife, skipping ", gene)
    next
  }

  base_dir <- file.path(PATHS$rates, "jackknife", gene)
  dir.create(base_dir, recursive = TRUE, showWarnings = FALSE)

  for (ds in DRAW_SIZES) {
    out_dir <- if (ds$label == "n1000") {
      base_dir
    } else {
      file.path(base_dir, ds$label)
    }
    res_out <- file.path(out_dir, "jackknife_results.csv")

    if (file.exists(res_out)) {
      message("  [", ds$label, "] Results already exist, skipping: ", res_out)
      next
    }

    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

    actual_size <- if (!is.na(ds$fixed)) {
      min(ds$fixed, n_avail)
    } else {
      round(n_avail * ds$frac)
    }

    message("  [", ds$label, "] Jackknife analysis (", N_SAMPLES, " × ",
            actual_size, " sequences)...")

    set.seed(MASTER_SEED)

    results_list    <- vector("list", N_SAMPLES)
    membership_list <- vector("list", N_SAMPLES)

    for (i in seq_len(N_SAMPLES)) {
      sampled_acc <- sample(df$accession, actual_size, replace = FALSE)

      sub_df   <- df |> filter(accession %in% sampled_acc)
      rownames(sub_df) <- sub_df$accession
      sub_tree <- unroot(keep.tip(tree, sampled_acc))

      fit <- tryCatch(
        phylolm(
          jk_formula,
          data  = sub_df,
          phy   = sub_tree,
          model = PARAMS$phylolm_method
        ),
        error = function(e) {
          message("  Sample ", i, " failed: ", conditionMessage(e))
          NULL
        }
      )

      if (is.null(fit)) next

      coef_mat <- summary(fit)$coefficients
      keep_terms <- intersect(PRED, rownames(coef_mat))
      if (length(keep_terms) == 0) next

      se_col <- intersect(
        c("StdErr", "Std.Error", "std_error"), colnames(coef_mat)
      )[1]
      t_col <- intersect(
        c("t.value", "t value", "z.value"), colnames(coef_mat)
      )[1]

      results_list[[i]] <- tibble(
        gene      = gene,
        sample_id = i,
        n         = actual_size,
        draw_size = ds$label,
        term      = keep_terms,
        estimate  = coef_mat[keep_terms, "Estimate"],
        std_error = coef_mat[keep_terms, se_col],
        t_value   = if (!is.na(t_col)) coef_mat[keep_terms, t_col] else NA_real_,
        p.value   = coef_mat[keep_terms, "p.value"],  # nolint: object_name_linter
        lambda    = fit$optpar
      )

      membership_list[[i]] <- tibble(
        sample_id = i,
        accession = sampled_acc
      )

      if (i %% 10 == 0) message("  Sample ", i, "/", N_SAMPLES, " done")
    }

    results_df    <- bind_rows(results_list)
    membership_df <- bind_rows(membership_list)

    write_csv(results_df,    res_out)
    write_csv(membership_df, file.path(out_dir, "jackknife_sequence_membership.csv"))

    message("  Written: ", out_dir)
    message("  Estimate range: [",
            round(min(results_df$estimate), 4), ", ",
            round(max(results_df$estimate), 4), "]")
  }
}

message("\nDone.")
