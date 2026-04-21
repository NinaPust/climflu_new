# scripts/00_subsample_test.R
#
# PURPOSE: Create a small stratified subsample of the pruned trees + matching
#          metadata for fast end-to-end testing of treetime and rates scripts.
#
#   Sampling is stratified by year so temporal signal is preserved.
#   Outputs mirror the real data structure under data/test/, which downstream
#   scripts pick up when passed the --test flag.
#
# USAGE:
#   Rscript scripts/00_subsample_test.R [--gene X] [--n 200] [--seed 42]
#
# INPUT:
#   data/trees/pruned/{gene}_pruned.nwk
#   data/metadata/{gene}_metadata.csv
#   data/climate/{gene}_climate.csv     (if present)
#
# OUTPUT:
#   data/test/trees/pruned/{gene}_pruned.nwk
#   data/test/metadata/{gene}_metadata.csv
#   data/test/climate/{gene}_climate.csv   (if source exists)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
  library(dplyr)
})

source(here("config.R"))

# ── Parse arguments ───────────────────────────────────────────────────────────
ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else ALL_GENES
n_target <- if ("--n" %in% args) as.integer(args[which(args == "--n") + 1]) else 200
seed     <- if ("--seed" %in% args) as.integer(args[which(args == "--seed") + 1]) else 42

set.seed(seed)
message("Subsampling ", n_target, " tips per gene (seed=", seed, ")")

# ── Create output dirs ────────────────────────────────────────────────────────
for (p in c(PATHS$test_trees_pruned, PATHS$test_trees_timed,
            PATHS$test_metadata, PATHS$test_climate)) {
  dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

for (gene in genes) {
  message("\n[", gene, "]")

  tree_in <- file.path(PATHS$trees_pruned, paste0(gene, "_pruned.nwk"))
  meta_in <- file.path(PATHS$metadata,     paste0(gene, "_metadata.csv"))

  if (!file.exists(tree_in)) { warning("  Pruned tree not found: ", tree_in); next }
  if (!file.exists(meta_in)) { warning("  Metadata not found: ",   meta_in); next }

  tree <- read.tree(tree_in)
  meta <- read_csv(meta_in, show_col_types = FALSE)

  # Restrict metadata to tips present in tree
  meta <- meta |> filter(accession %in% tree$tip.label)

  # Stratified sample by year — allocate proportionally, min 1 per year
  year_counts <- meta |> count(year) |> arrange(year)
  total_avail <- sum(year_counts$n)
  n_sample    <- min(n_target, total_avail)

  year_alloc <- year_counts |>
    mutate(alloc = pmax(1L, round(n * n_sample / total_avail))) |>
    mutate(alloc = pmin(alloc, n))

  # Trim allocation to hit n_sample exactly (remove from largest years first)
  diff <- sum(year_alloc$alloc) - n_sample
  if (diff > 0) {
    order_idx <- order(-year_alloc$alloc)
    for (i in seq_len(diff)) {
      idx <- order_idx[((i - 1) %% nrow(year_alloc)) + 1]
      if (year_alloc$alloc[idx] > 1) year_alloc$alloc[idx] <- year_alloc$alloc[idx] - 1L
    }
  }

  sampled_acc <- meta |>
    left_join(year_alloc |> select(year, alloc), by = "year") |>
    group_by(year) |>
    group_modify(~ slice_sample(.x, n = .x$alloc[1])) |>
    ungroup() |>
    pull(accession)

  # Prune tree and metadata to sample
  drop  <- setdiff(tree$tip.label, sampled_acc)
  stree <- drop.tip(tree, drop)
  smeta <- meta |> filter(accession %in% stree$tip.label)

  message("  ", Ntip(tree), " → ", Ntip(stree), " tips across ",
          n_distinct(smeta$year), " years (",
          min(smeta$year), "–", max(smeta$year), ")")

  write.tree(stree, file.path(PATHS$test_trees_pruned, paste0(gene, "_pruned.nwk")))
  write_csv(smeta, file.path(PATHS$test_metadata, paste0(gene, "_metadata.csv")))

  # Subset climate file if it exists
  climate_in <- file.path(PATHS$climate, paste0(gene, "_climate.csv"))
  if (file.exists(climate_in)) {
    climate <- read_csv(climate_in, show_col_types = FALSE)
    climate |>
      filter(accession %in% stree$tip.label) |>
      write_csv(file.path(PATHS$test_climate, paste0(gene, "_climate.csv")))
    message("  Climate subset written.")
  }
}

message("\nDone. Run downstream scripts with --test to use these subsets.")
