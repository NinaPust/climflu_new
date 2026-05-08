# scripts/02_phylogenetics/09_map_to_tree.R
#
# PURPOSE: Map assembled pairwise dN/dS matrices onto the MAPLE r2 pruned
#          tree topology using non-negative least squares (nnls.tree), giving
#          branch-specific dN and dS estimates for downstream rates analysis.
#
# INPUT:
#   data/rates/{gene}_dNdS_matrix.rds     — assembled dN/dS matrices
#   data/trees/pruned/{gene}_pruned.nwk   — MAPLE r2 pruned tree topology
#
# OUTPUT:
#   data/trees/{gene}_dN_tree.nwk
#   data/trees/{gene}_dS_tree.nwk
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(phangorn)
})

source(here("config.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

# Replace NAs in a symmetric distance matrix with row/column means.
impute_na <- function(mat) {
  na_mask <- is.na(mat)
  if (!any(na_mask)) return(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  for (i in which(rowSums(na_mask) > 0)) {
    bad <- which(na_mask[i, ])
    mat[i, bad] <- row_means[i]
    mat[bad, i] <- row_means[i]
  }
  diag(mat) <- 0
  mat
}

for (gene in genes) {
  message("\n[", gene, "] Mapping dN/dS matrices to tree...")

  rds_path  <- file.path(PATHS$rates, paste0(gene, "_dNdS_matrix.rds"))
  tree_path <- file.path(PATHS$trees_pruned, paste0(gene, "_pruned.nwk"))

  out_dn <- file.path(PATHS$trees, paste0(gene, "_dN_tree.nwk"))
  out_ds <- file.path(PATHS$trees, paste0(gene, "_dS_tree.nwk"))
  if (file.exists(out_dn) && file.exists(out_ds)) {
    message("  Both trees already exist, skipping.")
    next
  }

  if (!file.exists(rds_path)) {
    warning("Matrix RDS not found, skipping: ", rds_path); next
  }
  if (!file.exists(tree_path)) {
    warning("Pruned tree not found, skipping: ", tree_path); next
  }

  mats <- readRDS(rds_path)
  tree <- read.tree(tree_path)

  # Reconcile: keep only tips present in both matrix and tree
  common <- intersect(mats$accessions, tree$tip.label)
  n_drop_tree <- length(tree$tip.label) - length(common)
  n_drop_mat  <- length(mats$accessions) - length(common)
  message("  Tree tips: ", length(tree$tip.label),
          "  Matrix accessions: ", length(mats$accessions),
          "  Common: ", length(common))
  if (n_drop_tree > 0)
    message("  Dropping ", n_drop_tree, " tree tips not in matrix")
  if (n_drop_mat > 0)
    message("  Dropping ", n_drop_mat, " matrix accessions not in tree")

  tree   <- keep.tip(tree, common)
  dN_mat <- mats$dN[common, common]
  dS_mat <- mats$dS[common, common]

  # Impute missing pairwise distances with row means
  na_dN <- sum(is.na(dN_mat)) / 2
  na_dS <- sum(is.na(dS_mat)) / 2
  if (na_dN > 0) message("  Imputing ", na_dN, " missing dN pairs")
  if (na_dS > 0) message("  Imputing ", na_dS, " missing dS pairs")
  dN_mat <- impute_na(dN_mat)
  dS_mat <- impute_na(dS_mat)

  message("  Fitting dN tree (nnls.tree)...")
  t0 <- proc.time()
  dN_tree <- nnls.tree(as.dist(dN_mat), tree, rooted = FALSE)
  message("  Done (", round((proc.time() - t0)[["elapsed"]]), "s)")

  message("  Fitting dS tree (nnls.tree)...")
  t0 <- proc.time()
  dS_tree <- nnls.tree(as.dist(dS_mat), tree, rooted = FALSE)
  message("  Done (", round((proc.time() - t0)[["elapsed"]]), "s)")

  write.tree(dN_tree, out_dn)
  write.tree(dS_tree, out_ds)
  message("  dN tree: ", out_dn)
  message("  dS tree: ", out_ds)
}

message("\nDone.")
