# scripts/03_rates_analysis/02b_phylolm_dnds.R
#
# PURPOSE: Phylogenetic regression of per-branch dN rate, dS rate, and
#          branch omega (dN/dS) against climate predictors.
#          Uses fastme.bal terminal branch lengths divided by timetree
#          branch length as response variables.
#
# MODELS:
#   log(dN_rate | dS_rate | branch_omega) ~
#       scaled_temp + scaled_humidity + [+ scaled_pop_density]
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_phylolm_dnds_results.csv
#   data/rates/all_phylolm_dnds_results.csv   (FDR-corrected, all genes)
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
source(here("R", "phylolm_utils.R"))

args     <- commandArgs(trailingOnly = TRUE)
genes    <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES
use_test <- "--test" %in% args

if (use_test) message("*** Running in TEST mode (subsampled data) ***")

paths_timed <- if (use_test) PATHS$test_trees_timed else PATHS$trees_timed
paths_rates <- if (use_test) PATHS$test_rates       else PATHS$rates
paths_out   <- if (use_test) PATHS$test_rates       else PATHS$rates_phylolm

RESPONSES <- c("log_dN_rate", "log_dS_rate", "log_branch_omega")

PREDICTOR_TERMS <- c("scaled_temp", "scaled_humidity",
                     "scaled_delta_temp", "scaled_delta_humidity",
                     "scaled_sampling_density", "scaled_pop_density")

for (gene in genes) {
  message("\n[", gene, "] phylolm dN/dS analysis...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_out, paste0(gene, "_phylolm_dnds_results.csv"))

  if (!file.exists(rates_path)) {
    warning("Rates file not found: ", rates_path); next
  }
  if (!file.exists(tree_path)) {
    warning("Tree not found: ", tree_path); next
  }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  avail_responses <- intersect(RESPONSES, names(df))
  avail_responses <- avail_responses[
    sapply(avail_responses, function(r) sum(!is.na(df[[r]])) > 10)
  ]
  if (length(avail_responses) == 0) {
    message("  No dN/dS rate columns available, skipping.")
    next
  }

  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    message("  Dropping ", length(zero_tips), " zero-length tips")
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |>
    filter(accession %in% common) |>
    arrange(match(accession, tree$tip.label))
  rownames(df) <- df$accession

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

  message("  ", nrow(df), " taxa  |  responses: ",
          paste(avail_responses, collapse = ", "))

  all_results <- list()

  for (response in avail_responses) {
    sub_df   <- df |> filter(!is.na(.data[[response]]))
    sub_tree <- keep.tip(tree, sub_df$accession)

    if (nrow(sub_df) < 10) {
      message("  [", response, "] Too few observations, skipping.")
      next
    }

    year_term    <- if ("year" %in% names(sub_df)) "as.factor(year)" else NULL
    formula_full <- as.formula(
      paste(response, "~",
            paste(c(predictors, year_term), collapse = " + "))
    )

    result <- run_phylolm(
      formula  = formula_full,
      data     = sub_df,
      tree     = sub_tree,
      method   = PARAMS$phylolm_method,
      subset   = "full",
      response = response,
      gene     = gene
    )
    all_results[[paste(gene, response, sep = "_")]] <- result
  }

  results_df <- bind_rows(all_results)
  write_csv(results_df, out_path)
  message("  Written: ", out_path)
  print(results_df |> select(gene, response, term, estimate, p.value))
}

# ── FDR correction across all genes and responses ─────────────────────────────
message("\nApplying BH-FDR correction across all tests...")

all_res <- lapply(genes, function(g) {
  path <- file.path(paths_out, paste0(g, "_phylolm_dnds_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

if (nrow(all_res) > 0) {
  predictor_rows <- all_res$term %in% PREDICTOR_TERMS
  all_res$p.fdr  <- NA_real_
  all_res$p.fdr[predictor_rows] <- p.adjust(
    all_res$p.value[predictor_rows], method = "BH"
  )
  combined_path <- file.path(paths_out, "all_phylolm_dnds_results.csv")
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
}

message("\nDone.")
