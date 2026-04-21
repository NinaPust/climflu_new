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
                  if (has_pop) "scaled_pop_density",
                  if ("scaled_season" %in% names(df)) "scaled_season")

  message("  ", nrow(df), " taxa in analysis")

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

      formula_full <- as.formula(
        paste(response, "~", paste(predictors, collapse = " + "))
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

message("\nDone.")
