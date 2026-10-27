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
  library(parallel)
})

source(here("revisions", "config.R"))

out_dir <- file.path(PATHS$trees, "fastmebal")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

args     <- commandArgs(trailingOnly = TRUE)
genes    <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES
n_cores  <- min(length(genes) * 2L, detectCores() - 2L)

impute_na <- function(mat) {
  na_mask <- is.na(mat)
  if (!any(na_mask)) return(mat)
  row_means <- rowMeans(mat, na.rm = TRUE)
  for (i in which(rowSums(na_mask) > 0)) {
    bad         <- which(na_mask[i, ])
    mat[i, bad] <- row_means[i]
    mat[bad, i] <- row_means[i]
  }
  diag(mat) <- 0
  mat
}

build_one <- function(job) {
  gene   <- job$gene
  type   <- job$type   # "dN" or "dS"
  out_nwk <- job$out_nwk

  if (file.exists(out_nwk)) {
    message("[", gene, " ", type, "] Already exists, skipping.")
    return(invisible(NULL))
  }

  rds_path <- file.path(PATHS$rates, paste0(gene, "_dNdS_matrix.rds"))
  if (!file.exists(rds_path)) {
    warning("[", gene, " ", type, "] Matrix RDS not found: ", rds_path)
    return(invisible(NULL))
  }

  mats <- readRDS(rds_path)
  mat  <- mats[[type]][mats$accessions, mats$accessions]
  n_na <- sum(is.na(mat)) / 2
  if (n_na > 0)
    message("[", gene, " ", type, "] Imputing ", n_na, " missing pairs")
  mat <- impute_na(mat)

  message("[", gene, " ", type, "] Building tree (", length(mats$accessions), " taxa)...")
  t0   <- proc.time()
  tree <- fastme.bal(as.dist(mat))
  elapsed <- round((proc.time() - t0)[["elapsed"]])
  write.tree(tree, out_nwk)
  message("[", gene, " ", type, "] Done (", elapsed, "s) → ", out_nwk)
  invisible(NULL)
}

jobs <- lapply(genes, function(g) {
  list(
    list(gene = g, type = "dN",
         out_nwk = file.path(out_dir, paste0(g, "_dN_tree.nwk"))),
    list(gene = g, type = "dS",
         out_nwk = file.path(out_dir, paste0(g, "_dS_tree.nwk")))
  )
}) |> unlist(recursive = FALSE)

message("Launching ", length(jobs), " tree builds across ", n_cores, " cores...")
mclapply(jobs, build_one, mc.cores = n_cores)

message("\nDone.")
