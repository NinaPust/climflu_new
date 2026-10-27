# scripts/04_epitope_analysis/07_phylolm_epitope.R
#
# PURPOSE: Phylogenetic regression of epitope and non-epitope dN/dS rates
#          against climate predictors. Key comparison: does climate influence
#          selection differently at antigenic vs constrained sites?
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_phylolm_epitope_results.csv
#   data/rates/all_phylolm_epitope_results.csv
#
# FLAGS:
#   --gene X   run a single gene (default: H1 H3)

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "phylolm_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  c("H1", "H3")
}

predictor_terms <- c(
  "scaled_temp", "scaled_humidity",
  "scaled_delta_temp", "scaled_delta_humidity",
  "scaled_sampling_density", "scaled_pop_density"
)

# Response columns produced by 06_compute_epitope_rates.R
responses <- c(
  "log_epitope_dN_rate",    "log_epitope_dS_rate",    "log_epitope_omega",
  "log_nonepitope_dN_rate", "log_nonepitope_dS_rate", "log_nonepitope_omega",
  "log_ridge_dN_rate",      "log_ridge_dS_rate",      "log_ridge_omega",
  "log_nonridge_dN_rate",   "log_nonridge_dS_rate",   "log_nonridge_omega"
)

for (gene in genes) {
  message("\n[", gene, "] phylolm epitope analysis...")

  rates_path <- file.path(PATHS$rates,
                          paste0(gene, "_rates_epitope.csv"))
  tree_path  <- file.path(PATHS$trees_timed,
                          paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(PATHS$rates,
                          paste0(gene, "_phylolm_epitope_results.csv"))

  if (!file.exists(rates_path)) {
    warning("Rates file not found: ", rates_path); next
  }
  if (!file.exists(tree_path)) {
    warning("Timetree not found: ", tree_path); next
  }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  avail_responses <- intersect(responses, names(df))
  avail_responses <- avail_responses[
    sapply(avail_responses,
           function(r) sum(!is.na(df[[r]])) > 10)
  ]
  if (length(avail_responses) == 0) {
    message("  No epitope rate columns available, skipping."); next
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

  has_pop    <- "scaled_pop_density" %in% names(df)
  predictors <- c(
    "scaled_temp", "scaled_humidity",
    if ("scaled_delta_temp"       %in% names(df)) "scaled_delta_temp",
    if ("scaled_delta_humidity"   %in% names(df)) "scaled_delta_humidity",
    if ("scaled_sampling_density" %in% names(df)) "scaled_sampling_density",
    if (has_pop) "scaled_pop_density"
  )

  message("  ", nrow(df), " taxa  |  responses: ",
          paste(avail_responses, collapse = ", "))

  all_results <- list()

  for (response in avail_responses) {
    sub_df   <- df |> filter(!is.na(.data[[response]]))
    sub_tree <- keep.tip(tree, sub_df$accession)
    if (nrow(sub_df) < 10) next

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
  print(results_df |>
          filter(term %in% predictor_terms) |>
          select(gene, response, term, estimate, p.value))
}

# ── BH-FDR correction across all genes ────────────────────────────────────────
message("\nApplying BH-FDR correction...")

all_res <- lapply(genes, function(g) {
  path <- file.path(PATHS$rates, paste0(g, "_phylolm_epitope_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

if (nrow(all_res) > 0) {
  pred_rows      <- all_res$term %in% predictor_terms
  all_res$p.fdr  <- NA_real_
  all_res$p.fdr[pred_rows] <- p.adjust(
    all_res$p.value[pred_rows], method = "BH"
  )
  combined_path <- file.path(PATHS$rates, "all_phylolm_epitope_results.csv")
  write_csv(all_res, combined_path)
  message("  Written: ", combined_path)

  sig <- all_res |>
    filter(term %in% predictor_terms, p.fdr < 0.05) |>
    select(gene, response, term, estimate, p.value, p.fdr) |>
    mutate(across(c(estimate), \(x) round(x, 4)),
           across(c(p.value, p.fdr), \(x) signif(x, 3))) |>
    arrange(p.fdr)
  message("\nSignificant after FDR (p.fdr < 0.05):")
  print(sig)
}

message("\nDone.")
