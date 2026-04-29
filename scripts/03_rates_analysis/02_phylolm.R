# scripts/03_rates_analysis/02_phylolm.R
#
# PURPOSE: Run phylogenetic comparative regression (phylolm) to test whether
#          temperature and humidity predict evolutionary rates and dN/dS,
#          accounting for phylogenetic non-independence.
#
# MODELS:
#   log(rate|omega) ~ scaled_temp + scaled_humidity + scaled_pop_density
#   Tips with NA pop_density are dropped before fitting.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_phylolm_results.csv  — coefficient table (all models)
#
# FLAGS:
#   --gene X   run a single gene segment
#   --test     use subsampled data from data/test/

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
  library(MuMIn)
  library(tidyr)
})

source(here("config.R"))
source(here("R", "phylolm_utils.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args     <- commandArgs(trailingOnly = TRUE)
genes    <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}
use_test <- "--test" %in% args

if (use_test) message("*** Running in TEST mode (subsampled data) ***")

paths_timed <- if (use_test) PATHS$test_trees_timed else PATHS$trees_timed
paths_rates <- if (use_test) PATHS$test_rates       else PATHS$rates

dir.create(paths_rates, recursive = TRUE, showWarnings = FALSE)

for (gene in genes) {
  message("\n[", gene, "] phylolm analysis...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_rates, paste0(gene, "_phylolm_results.csv"))

  if (!file.exists(rates_path)) {
    warning("Rates file not found: ", rates_path)
    next
  }
  if (!file.exists(tree_path)) {
    warning("Tree not found: ", tree_path)
    next
  }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  # Align tips and data
  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  # Drop tips with zero-length terminal branches — they cause a singular
  # covariance matrix in phylolm ("sister external edges have length 0").
  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    message("  Dropping ", length(zero_tips),
            " tips with zero-length branches (phylolm singularity)")
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |>
    filter(accession %in% common) |>
    arrange(match(accession, tree$tip.label))
  rownames(df) <- df$accession

  # Drop tips with NA pop_density
  has_pop <- "scaled_pop_density" %in% names(df)
  if (has_pop) {
    n_before <- nrow(df)
    df   <- df |> filter(!is.na(scaled_pop_density))
    tree <- keep.tip(tree, df$accession)
    if (nrow(df) < n_before)
      message("  Dropped ", n_before - nrow(df), " tips with NA pop_density")
  }

  predictors <- c("scaled_temp", "scaled_humidity",
                  if ("scaled_delta_temp"       %in% names(df)) "scaled_delta_temp",
                  if ("scaled_delta_humidity"   %in% names(df)) "scaled_delta_humidity",
                  if ("scaled_sampling_density" %in% names(df)) "scaled_sampling_density",
                  if (has_pop) "scaled_pop_density")
                  # if ("scaled_season" %in% names(df)) "scaled_season"

  message("  ", nrow(df), " taxa in analysis")

  # ── VIF (computed from OLS — phylogenetic structure doesn't affect collinearity)
  if (length(predictors) > 1) {
    vif_vals <- sapply(predictors, function(p) {
      others <- setdiff(predictors, p)
      f  <- as.formula(paste(p, "~", paste(others, collapse = " + ")))
      r2 <- summary(lm(f, data = df))$r.squared
      round(1 / (1 - r2), 2)
    })
    message("  VIF: ", paste(names(vif_vals), vif_vals, sep = " = ",
                              collapse = "  |  "))
  }

  all_results <- list()

  for (cfg in list(list(subset_name = "full", data = df))) {
    sub_df   <- cfg$data
    sub_name <- cfg$subset_name

    if (nrow(sub_df) < 10) {
      message("  [", sub_name, "] Too few observations (",
              nrow(sub_df), "), skipping.")
      next
    }
    sub_tree <- keep.tip(tree, sub_df$accession)

    for (response in c("log_rate", "log_omega")) {
      if (!response %in% names(sub_df) || all(is.na(sub_df[[response]]))) next

      # Year factor controls for secular trends (not included in predictors so
      # it is excluded from FDR correction and forest plots).
      year_term <- if ("year" %in% names(sub_df)) "as.factor(year)" else NULL
      formula_full <- as.formula(
        paste(response, "~",
              paste(c(predictors, year_term), collapse = " + "))
      )

      result <- run_phylolm(
        formula  = formula_full,
        data     = sub_df,
        tree     = sub_tree,
        method   = PARAMS$phylolm_method,
        subset   = sub_name,
        response = response,
        gene     = gene
      )
      all_results[[paste(gene, sub_name, response, sep = "_")]] <- result
    }
  }

  results_df <- bind_rows(all_results)
  write_csv(results_df, out_path)
  message("  Written: ", out_path)
  print(results_df |> select(gene, subset, response, term, estimate, p.value))
}

# ── FDR correction across all genes and responses ─────────────────────────────
# Apply Benjamini-Hochberg correction over all predictor terms (not intercepts)
# across all genes and response variables. Per-gene CSVs are left unchanged;
# corrected results go to a single combined file.
message("\nApplying BH-FDR correction across all tests...")

PREDICTOR_TERMS <- c("scaled_temp", "scaled_humidity",
                     "scaled_delta_temp", "scaled_delta_humidity",
                     "scaled_sampling_density", "scaled_pop_density")
                     # "scaled_season"

all_res <- lapply(genes, function(g) {
  path <- file.path(paths_rates, paste0(g, "_phylolm_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

predictor_rows <- all_res$term %in% PREDICTOR_TERMS
all_res$p.fdr <- NA_real_
all_res$p.fdr[predictor_rows] <- p.adjust(
  all_res$p.value[predictor_rows], method = "BH"
)

combined_path <- file.path(paths_rates, "all_phylolm_results.csv")
write_csv(all_res, combined_path)
message("  Written: ", combined_path)

message("\nSignificant after FDR (p.fdr < 0.05):")
sig <- all_res |>
  filter(term %in% PREDICTOR_TERMS, p.fdr < 0.05) |>
  select(gene, response, term, estimate, p.value, p.fdr) |>
  mutate(across(c(estimate), \(x) round(x, 4)),
         across(c(p.value, p.fdr), \(x) signif(x, 3))) |>
  arrange(p.fdr)
print(sig)

message("\nCorrelation matrix of predictors (first gene with full data):")
ref_gene <- genes[1]
ref_path <- file.path(paths_rates, paste0(ref_gene, "_rates.csv"))
if (file.exists(ref_path)) {
  ref_df <- read_csv(ref_path, show_col_types = FALSE)
  avail  <- intersect(PREDICTOR_TERMS, names(ref_df))
  if (length(avail) > 1) {
    cormat <- ref_df |>
      select(all_of(avail)) |>
      drop_na() |>
      cor() |>
      round(2)
    print(cormat)
  }
}

message("\nDone.")
