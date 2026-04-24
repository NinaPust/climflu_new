# scripts/02_phylogenetics/08_pagels_lambda.R
#
# PURPOSE: Estimate Pagel's lambda (phylogenetic signal in log_rate) on
#          MAPLE trees from each rebuild round to assess whether outlier
#          removal improved phylogenetic signal.
#
# APPROACH:
#   For each gene × round, 1000 tips are stratified-subsampled (by year),
#   the tree is pruned to that subset, and phylosig() fits lambda via ML.
#   The same seed and subsample size as 04_inla_phylo.R are used so results
#   are directly comparable to the INLA analysis.
#
#   Rounds compared:
#     r1 — pruned r1 MAPLE tree  (data/trees/pruned/{gene}_pruned.nwk)
#     r2 — raw r2 MAPLE tree     (data/trees/raw/maple/{gene}_tree_r2.nwk)
#     r3 — raw r3 MAPLE tree     (data/trees/raw/maple/{gene}_tree_r3.nwk)
#          (only for genes where r3 was built)
#
# INPUT:
#   data/trees/pruned/{gene}_pruned.nwk
#   data/trees/raw/maple/{gene}_tree_r2.nwk
#   data/trees/raw/maple/{gene}_tree_r3.nwk   (if present)
#   data/rates/{gene}_rates.csv
#
# OUTPUT:
#   data/trees/pagels_lambda.csv
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(phytools)
  library(dplyr)
  library(readr)
})

source(here("config.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args      <- commandArgs(trailingOnly = TRUE)
genes     <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}

N_SUB <- 1000

stratified_subsample <- function(df, n, seed = 42) {
  set.seed(seed)
  if (nrow(df) <= n) return(df)
  total_n    <- nrow(df)
  year_counts <- df |>
    count(.data$year, name = "grp_n") |>
    mutate(.k = pmax(1L, as.integer(round(.data$grp_n / total_n * n))))
  sampled <- lapply(split(df, df$year), function(sub) {
    k <- year_counts$.k[year_counts$year == sub$year[1]]
    sub[sample(nrow(sub), min(k, nrow(sub))), ]
  })
  result <- do.call(rbind, sampled)
  result[sample(nrow(result), min(n, nrow(result))), ]
}

out_path <- file.path(PATHS$trees, "pagels_lambda.csv")

all_results <- list()

for (gene in genes) {

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))
  if (!file.exists(rates_path)) {
    warning("Rates file not found, skipping: ", rates_path)
    next
  }
  df <- read_csv(rates_path, show_col_types = FALSE) |>
    filter(is.finite(log_rate))
  if ("scaled_pop_density" %in% names(df))
    df <- df |> filter(!is.na(scaled_pop_density))

  round_paths <- list(
    r1 = file.path(PATHS$trees_pruned,
                   paste0(gene, "_pruned_r1.nwk")),
    r2 = file.path(PATHS$trees_pruned,
                   paste0(gene, "_pruned.nwk"))
  )

  for (round in names(round_paths)) {
    tree_path <- round_paths[[round]]

    if (!file.exists(tree_path)) {
      message("[", gene, " ", round, "] No tree found, skipping.")
      next
    }

    message("\n[", gene, " ", round, "] Loading tree...")
    tree <- read.tree(tree_path)

    # Drop zero-length tips (duplicates) — mirrors 04_inla_phylo.R
    tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
    zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
    if (length(zero_tips) > 0) tree <- drop.tip(tree, zero_tips)

    common <- intersect(tree$tip.label, df$accession)
    tree   <- keep.tip(tree, common)
    sub_df <- df |> filter(accession %in% tree$tip.label)

    message("  ", nrow(sub_df), " taxa available; subsampling to ",
            min(N_SUB, nrow(sub_df)))

    sub_df   <- stratified_subsample(sub_df, N_SUB)
    sub_tree <- keep.tip(tree, sub_df$accession)

    x <- setNames(sub_df$log_rate, sub_df$accession)

    message("  Fitting Pagel's lambda...")
    res <- tryCatch(
      phylosig(sub_tree, x, method = "lambda", test = TRUE),
      error = function(e) {
        warning("phylosig failed for ", gene, " ", round, ": ", e$message)
        NULL
      }
    )
    if (is.null(res)) next

    row <- tibble(
      gene   = gene,
      round  = round,
      n      = length(x),
      lambda = res$lambda,
      logL   = res$logL,
      logL0  = res$logL0,
      p      = res$P
    )
    message(sprintf("  lambda=%.4f  logL=%.2f  logL0=%.2f  p=%.4f",
                    res$lambda, res$logL, res$logL0, res$P))

    all_results[[paste(gene, round)]] <- row
  }
}

if (length(all_results) > 0) {
  combined <- bind_rows(all_results)
  write_csv(combined, out_path)
  message("\nSaved: ", out_path)
  message("\n── Results ─────────────────────────────────────────────────")
  print(combined, n = 50)
} else {
  message("No results to save.")
}

message("\nDone.")
