# scripts/03_rates_analysis/02_phylolm.R
#
# PURPOSE: Run phylogenetic comparative regression (phylolm) to test whether
#          temperature and humidity predict evolutionary rates and dN/dS,
#          accounting for phylogenetic non-independence.
#
# MODELS:
#   log(rate)  ~ scaled_temp + scaled_humidity  [+ pop_density if available]
#   log(omega) ~ scaled_temp + scaled_humidity
#   Both models run on full data, then urban-only, rural-only subsets.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_phylolm_results.csv  — coefficient table (all models)

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

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] phylolm analysis...")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(PATHS$trees, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(PATHS$rates, paste0(gene, "_phylolm_results.csv"))

  if (!file.exists(rates_path)) { warning("Rates file not found: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Tree not found: ", tree_path);         next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  # Align tips and data
  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)
  df     <- df |> filter(accession %in% common) |>
    arrange(match(accession, tree$tip.label))
  rownames(df) <- df$accession

  message("  ", length(common), " taxa in analysis")

  # ── Define subsets ────────────────────────────────────────────────────────
  # Urban/rural split is approximate — refine thresholds as needed
  subsets <- list(
    full  = df,
    urban = df |> filter(!is.na(pop_density) & pop_density > 500),
    rural = df |> filter(!is.na(pop_density) & pop_density <= 500)
  )

  all_results <- list()

  for (subset_name in names(subsets)) {
    sub_df <- subsets[[subset_name]]
    if (nrow(sub_df) < 10) {
      message("  [", subset_name, "] Too few observations (", nrow(sub_df), "), skipping.")
      next
    }
    sub_tree <- keep.tip(tree, sub_df$accession)

    for (response in c("log_rate", "log_omega")) {
      if (!response %in% names(sub_df) || all(is.na(sub_df[[response]]))) next

      predictors <- intersect(c("scaled_temp", "scaled_humidity", "scaled_pop_density"),
                               names(sub_df))

      formula_full <- as.formula(paste(response, "~", paste(predictors, collapse = " + ")))

      result <- run_phylolm(
        formula  = formula_full,
        data     = sub_df,
        tree     = sub_tree,
        method   = PARAMS$phylolm_method,
        subset   = subset_name,
        response = response,
        gene     = gene
      )
      all_results[[paste(gene, subset_name, response, sep = "_")]] <- result
    }
  }

  results_df <- bind_rows(all_results)
  write_csv(results_df, out_path)
  message("  Written: ", out_path)
  print(results_df |> select(gene, subset, response, term, estimate, p.value))
}

message("\nDone.")
