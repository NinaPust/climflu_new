# scripts/03_rates_analysis/07_phylolm_jackknife.R
#
# PURPOSE: Block jackknife subsampling for phylolm temperature analysis.
#          Runs 100 random samples at four draw sizes (1000 fixed, 15%, 30%,
#          60% of dataset), runs phylolm (log_rate ~ scaled_temp +
#          as.factor(year)) on each, and saves per-gene results for downstream
#          robustness plots.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
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

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))

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

      se_col <- intersect(
        c("StdErr", "Std.Error", "std_error"), colnames(coef_mat)
      )[1]

      results_list[[i]] <- tibble(
        gene      = gene,
        sample_id = i,
        n         = actual_size,
        draw_size = ds$label,
        estimate  = coef_mat["scaled_temp", "Estimate"],
        std_error = coef_mat["scaled_temp", se_col],
        t_value   = coef_mat["scaled_temp", "t.value"],
        p.value   = coef_mat["scaled_temp", "p.value"],  # nolint: object_name_linter
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
