# scripts/02_phylogenetics/01_maple_tree.R
#
# PURPOSE: Infer an initial phylogenetic tree with MAPLE, then prune outlier
#          taxa in two rounds based on branch length thresholds.
#
# INPUT:
#   data/alignments/{gene}_aligned.fasta
#
# OUTPUT:
#   data/trees/{gene}_tree.nwk                — pruned tree (Newick)
#   data/metadata/{gene}_accessions_final.txt — accessions retained after pruning

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("config.R"))
source(here("R", "tree_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] MAPLE tree inference + pruning...")

  fasta_in  <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  tree_out  <- file.path(PATHS$trees,      paste0(gene, "_tree.nwk"))
  acc_out   <- file.path(PATHS$metadata,   paste0(gene, "_accessions_final.txt"))
  maple_tmp <- file.path(PATHS$trees,      paste0(gene, "_maple_tmp"))

  if (!file.exists(fasta_in)) {
    warning("Alignment not found, skipping: ", fasta_in)
    next
  }

  # ── Run MAPLE ──────────────────────────────────────────────────────────────
  dir.create(maple_tmp, showWarnings = FALSE, recursive = TRUE)
  maple_cmd <- paste(TOOLS$maple, "--input", fasta_in, "--output", maple_tmp)
  ret <- system(maple_cmd)
  if (ret != 0) stop("MAPLE failed for gene ", gene)

  # MAPLE outputs a Newick file — locate it
  maple_tree_file <- list.files(maple_tmp, pattern = "\\.nwk$|\\.tree$", full.names = TRUE)[1]
  if (is.na(maple_tree_file)) stop("MAPLE output tree not found in: ", maple_tmp)

  tree <- read.tree(maple_tree_file)
  message("  Initial tree: ", Ntip(tree), " tips")

  # ── Round 1: remove tips with very short branches ─────────────────────────
  tree <- prune_outliers(tree,
                         min_bl = PARAMS$prune_min_bl,
                         max_bl = Inf,
                         round  = 1)

  # ── Round 2: remove tips with very long branches ──────────────────────────
  tree <- prune_outliers(tree,
                         min_bl = 0,
                         max_bl = PARAMS$prune_max_bl,
                         round  = 2)

  message("  Final tree: ", Ntip(tree), " tips")

  write.tree(tree, tree_out)
  writeLines(tree$tip.label, acc_out)
  message("  Tree:        ", tree_out)
  message("  Accessions:  ", acc_out)
}

message("\nDone.")
