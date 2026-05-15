# scripts/03_rates_analysis/07_phylolm_jackknife.R
#
# PURPOSE: Block jackknife subsampling for phylolm temperature analysis.
#          Takes 100 random samples of 1000 sequences, runs phylolm
#          (log_rate ~ scaled_temp + as.factor(year)) on each, and saves
#          per-gene results for downstream robustness plots.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/jackknife/{gene}/jackknife_results.csv
#   data/rates/jackknife/{gene}/jackknife_sequence_membership.csv
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
SAMPLE_SIZE <- 1000L
MASTER_SEED <- 42L

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  GENES
}

for (gene in genes) {
  message("\n[", gene, "] Jackknife analysis (", N_SAMPLES, " × ",
          SAMPLE_SIZE, " sequences)...")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
  out_dir    <- file.path(PATHS$rates, "jackknife", gene)

  if (!file.exists(rates_path)) { warning("Rates file not found: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Tree not found: ", tree_path);         next }

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

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
    filter(!is.na(scaled_temp), !is.na(log_rate)) |>
    arrange(match(accession, tree$tip.label))
  tree <- keep.tip(tree, df$accession)
  rownames(df) <- df$accession

  n_avail <- nrow(df)
  message("  ", n_avail, " taxa available for subsampling")

  if (n_avail < 10) {
    warning("  Too few taxa for jackknife, skipping ", gene)
    next
  }

  actual_size <- min(SAMPLE_SIZE, n_avail)
  if (actual_size < SAMPLE_SIZE) {
    warning("  Only ", n_avail, " sequences available; sampling all of them")
  }

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
        log_rate ~ scaled_temp + as.factor(year),
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
    if (!"scaled_temp" %in% rownames(coef_mat)) next

    se_col  <- intersect(c("StdErr", "Std.Error", "std_error"), colnames(coef_mat))[1]

    results_list[[i]] <- tibble(
      gene      = gene,
      sample_id = i,
      n         = actual_size,
      estimate  = coef_mat["scaled_temp", "Estimate"],
      std_error = coef_mat["scaled_temp", se_col],
      t_value   = coef_mat["scaled_temp", "t.value"],
      p.value   = coef_mat["scaled_temp", "p.value"],
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

  write_csv(results_df,    file.path(out_dir, "jackknife_results.csv"))
  write_csv(membership_df, file.path(out_dir, "jackknife_sequence_membership.csv"))

  message("  Written: ", out_dir)
  message("  Estimate range: [",
          round(min(results_df$estimate), 4), ", ",
          round(max(results_df$estimate), 4), "]")
}

message("\nDone.")
