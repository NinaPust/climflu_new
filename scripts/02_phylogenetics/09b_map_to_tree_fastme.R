# scripts/02_phylogenetics/09b_map_to_tree_fastme.R
#
# PURPOSE: Build dN and dS trees from assembled pairwise matrices using
#          fastme.bal (balanced minimum evolution). Faster alternative to
#          nnls.tree (09_map_to_tree.R) for large datasets (10k+ tips).
#          Terminal branch lengths are used downstream to compute per-lineage
#          dN/dS rates; topology is matched by tip label, not tree structure.
#
# INPUT:
#   data/rates/{gene}_dNdS_matrix.rds
#
# OUTPUT:
#   data/trees/fastmebal/{gene}_dN_tree.nwk
#   data/trees/fastmebal/{gene}_dS_tree.nwk
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("config.R"))

out_dir <- file.path(PATHS$trees, "fastmebal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

# Replace NAs with row means (symmetric imputation).
impute_na <- function(mat) {
  na_mask <- is.na(mat)
  if (!any(na_mask)) return(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  for (i in which(rowSums(na_mask) > 0)) {
    bad        <- which(na_mask[i, ])
    mat[i, bad] <- row_means[i]
    mat[bad, i] <- row_means[i]
  }
  diag(mat) <- 0
  mat
}

for (gene in genes) {
  message("\n[", gene, "] Building dN/dS trees with fastme.bal...")

  rds_path <- file.path(PATHS$rates, paste0(gene, "_dNdS_matrix.rds"))
  if (!file.exists(rds_path)) {
    warning("Matrix RDS not found, skipping: ", rds_path); next
  }

  mats <- readRDS(rds_path)
  message("  ", length(mats$accessions), " accessions")

  dN_mat <- mats$dN[mats$accessions, mats$accessions]
  dS_mat <- mats$dS[mats$accessions, mats$accessions]

  na_dN <- sum(is.na(dN_mat)) / 2
  na_dS <- sum(is.na(dS_mat)) / 2
  if (na_dN > 0) message("  Imputing ", na_dN, " missing dN pairs")
  if (na_dS > 0) message("  Imputing ", na_dS, " missing dS pairs")
  dN_mat <- impute_na(dN_mat)
  dS_mat <- impute_na(dS_mat)

  message("  Building dN tree...")
  t0 <- proc.time()
  dN_tree <- fastme.bal(as.dist(dN_mat))
  message("  Done (", round((proc.time() - t0)[["elapsed"]]), "s)")

  message("  Building dS tree...")
  t0 <- proc.time()
  dS_tree <- fastme.bal(as.dist(dS_mat))
  message("  Done (", round((proc.time() - t0)[["elapsed"]]), "s)")

  out_dn <- file.path(out_dir, paste0(gene, "_dN_tree.nwk"))
  out_ds <- file.path(out_dir, paste0(gene, "_dS_tree.nwk"))
  write.tree(dN_tree, out_dn)
  write.tree(dS_tree, out_ds)
  message("  dN tree: ", out_dn)
  message("  dS tree: ", out_ds)
}

message("\nDone.")
