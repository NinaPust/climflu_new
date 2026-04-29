# scripts/03_rates_analysis/06_phylolm_country.R
#
# PURPOSE: Sensitivity analysis — restrict to countries with dense year-round
#          surveillance, then add country fixed effects on top of the sampling
#          density and year FE controls already used in 02_phylolm.R.
#          This is robustness check 3 from the confound investigation.
#
# COUNTRIES: USA, United Kingdom, Australia, China, Japan.
#            Edit DENSE_COUNTRIES below to adjust.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_phylolm_country_results.csv
#   data/rates/all_phylolm_country_results.csv
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
  library(tidyr)
})

source(here("config.R"))
source(here("R", "phylolm_utils.R"))

DENSE_COUNTRIES <- c("USA", "United Kingdom", "Australia", "China", "Japan")

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
                     "scaled_delta_temp", "scaled_delta_humidity",
                     "scaled_sampling_density", "scaled_pop_density")
                     # "scaled_season"

for (gene in genes) {
  message("\n[", gene, "] phylolm (dense-country subset)...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_rates,
                          paste0(gene, "_phylolm_country_results.csv"))

  if (!file.exists(rates_path)) { warning("Missing: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Missing: ", tree_path);  next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  # ── Restrict to dense-surveillance countries ────────────────────────────
  n_all <- nrow(df)
  df    <- df |> filter(country %in% DENSE_COUNTRIES)
  message("  Country filter: ", n_all, " → ", nrow(df), " rows (",
          paste(sort(unique(df$country)), collapse = ", "), ")")

  if (nrow(df) < 10) {
    warning("Too few rows after country filter for ", gene, ", skipping.")
    next
  }

  # ── Standard tree/data alignment ────────────────────────────────────────
  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |>
    filter(accession %in% common) |>
    arrange(match(accession, tree$tip.label))

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

  # Drop country and year levels with fewer than 2 observations so factor
  # contrasts don't cause rank-deficiency.
  df <- df |>
    group_by(country) |> filter(n() >= 2) |> ungroup() |>
    group_by(year)    |> filter(n() >= 2) |> ungroup()
  tree <- keep.tip(tree, df$accession)

  message("  ", nrow(df), " taxa after level pruning")

  # Set rownames after all filtering — phylolm matches data to tips by rowname.
  # Must be done last since dplyr pipes silently drop rownames on tibbles.
  df <- as.data.frame(df)
  rownames(df) <- df$accession

  all_results <- list()

  for (response in c("log_rate", "log_omega")) {
    if (!response %in% names(df) || all(is.na(df[[response]]))) next

    formula_str <- paste(
      response, "~",
      paste(c(predictors,
              "as.factor(year)",
              "as.factor(country)"),
            collapse = " + ")
    )

    result <- run_phylolm(
      formula  = as.formula(formula_str),
      data     = df,
      tree     = tree,
      method   = PARAMS$phylolm_method,
      subset   = "dense_country",
      response = response,
      gene     = gene
    )
    all_results[[paste(gene, response, sep = "_")]] <- result
  }

  results_df <- bind_rows(all_results)
  write_csv(results_df, out_path)
  message("  Written: ", out_path)
  results_df |>
    filter(term %in% PREDICTOR_TERMS) |>
    select(gene, response, term, estimate, p.value) |>
    print()
}

# ── FDR correction ────────────────────────────────────────────────────────────
message("\nApplying BH-FDR correction...")

all_res <- lapply(genes, function(g) {
  path <- file.path(paths_rates, paste0(g, "_phylolm_country_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

predictor_rows    <- all_res$term %in% PREDICTOR_TERMS
all_res$p.fdr     <- NA_real_
all_res$p.fdr[predictor_rows] <- p.adjust(
  all_res$p.value[predictor_rows], method = "BH"
)

combined_path <- file.path(paths_rates, "all_phylolm_country_results.csv")
write_csv(all_res, combined_path)
message("  Written: ", combined_path)

message("\nSignificant after FDR (p.fdr < 0.05):")
all_res |>
  filter(term %in% PREDICTOR_TERMS, p.fdr < 0.05) |>
  select(gene, response, term, estimate, p.value, p.fdr) |>
  mutate(across(c(estimate), \(x) round(x, 4)),
         across(c(p.value, p.fdr), \(x) signif(x, 3))) |>
  arrange(p.fdr) |>
  print()

message("\nDone.")
