# scripts/03_rates_analysis/02c_phylolm_pandemic.R
#
# PURPOSE: Repeat the main phylolm (log_rate ~ climate + year) for H1 only,
#          split into pandemic (2009–2010) and post-pandemic (2011+) subsets.
#          as.factor(year) in the main model already absorbs secular trends,
#          but this supplementary analysis shows the climate signal holds
#          within each era independently.
#
# INPUT:
#   data/rates/H1_rates.csv
#   data/trees/timed/H1_timetree.nwk
#
# OUTPUT:
#   data/rates/phylolm/H1_phylolm_pandemic_results.csv

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
})

source(here("config.R"))
source(here("R", "phylolm_utils.R"))

gene      <- "H1"
rates_path <- file.path(PATHS$rates,       paste0(gene, "_rates.csv"))
tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
out_path   <- file.path(PATHS$rates_phylolm, "H1_phylolm_pandemic_results.csv")

df   <- read_csv(rates_path, show_col_types = FALSE)
tree <- read.tree(tree_path)

common <- intersect(tree$tip.label, df$accession)
tree   <- keep.tip(tree, common)

tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
if (length(zero_tips) > 0) {
  message("Dropping ", length(zero_tips), " zero-length tips")
  tree   <- drop.tip(tree, zero_tips)
}

df <- df |>
  filter(accession %in% tree$tip.label, !is.na(log_rate)) |>
  arrange(match(accession, tree$tip.label))
rownames(df) <- df$accession

predictors <- c(
  "scaled_temp", "scaled_humidity",
  if ("scaled_delta_temp"       %in% names(df)) "scaled_delta_temp",
  if ("scaled_delta_humidity"   %in% names(df)) "scaled_delta_humidity",
  if ("scaled_sampling_density" %in% names(df)) "scaled_sampling_density"
)

subsets <- list(
  pandemic      = df |> filter(year %in% c(2009, 2010)),
  post_pandemic = df |> filter(year >= 2011)
)

all_results <- list()

for (subset_name in names(subsets)) {
  sub_df <- subsets[[subset_name]]
  message("[H1 / ", subset_name, "] n = ", nrow(sub_df))

  if (nrow(sub_df) < 10) {
    message("  Too few observations, skipping.")
    next
  }

  sub_tree <- keep.tip(tree, sub_df$accession)

  # Only include year factor if more than one unique year in subset
  year_term <- if (length(unique(sub_df$year)) > 1) "as.factor(year)" else NULL
  formula_full <- as.formula(
    paste("log_rate ~", paste(c(predictors, year_term), collapse = " + "))
  )

  result <- run_phylolm(
    formula  = formula_full,
    data     = sub_df,
    tree     = sub_tree,
    method   = PARAMS$phylolm_method,
    subset   = subset_name,
    response = "log_rate",
    gene     = gene
  )
  all_results[[subset_name]] <- result
}

results_df <- bind_rows(all_results)

# FDR correction within this analysis
PREDICTOR_TERMS <- c("scaled_temp", "scaled_humidity", "scaled_delta_temp",
                     "scaled_delta_humidity", "scaled_sampling_density")
pred_rows <- results_df$term %in% PREDICTOR_TERMS
results_df$p.fdr <- NA_real_
results_df$p.fdr[pred_rows] <- p.adjust(results_df$p.value[pred_rows], method = "BH")

write_csv(results_df, out_path)
message("Written: ", out_path)

print(results_df |>
        filter(term %in% PREDICTOR_TERMS) |>
        select(subset, term, estimate, std_error, p.value, p.fdr))