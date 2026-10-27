# scripts/04_epitope_analysis/05b_map_to_tree_fastme.R
#
# PURPOSE: Build dN and dS trees from epitope and non-epitope pairwise
#          matrices using fastme.bal (balanced minimum evolution).
#          Mirrors 02_phylogenetics/09b_map_to_tree_fastme.R but operates
#          on data/epitope/ matrices. All 8 trees are run in parallel.
#          Terminal branch lengths from these trees are used in
#          06_compute_epitope_rates.R to derive per-lineage dN/dS rates.
#
# INPUT:   data/epitope/{gene}_{type}_dNdS_matrix.rds
# OUTPUT:  data/epitope/trees/{gene}_{type}_{dN|dS}_tree.nwk
#
# FLAGS:
#   --gene X   run a single gene (default: H1 H3)
#   --type X   epitope|nonepitope (default: both)

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("revisions", "config.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  c("H1", "H3")
}
types <- if ("--type" %in% args) {
  args[which(args == "--type") + 1]
} else {
  c("epitope", "nonepitope", "ridge", "nonridge")
}

epi_dir  <- here("revisions", "data", "epitope")
tree_dir <- file.path(epi_dir, "trees")
dir.create(tree_dir, showWarnings = FALSE, recursive = TRUE)

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
  gene     <- job$gene
  type     <- job$type
  mat_type <- job$mat_type   # "dN" or "dS"
  out_nwk  <- job$out_nwk

  if (file.exists(out_nwk)) {
    message("[", gene, " ", type, " ", mat_type, "] Already exists, skipping.")
    return(invisible(NULL))
  }

  rds_path <- file.path(epi_dir,
                        paste0(gene, "_", type, "_dNdS_matrix.rds"))
  if (!file.exists(rds_path)) {
    warning("[", gene, " ", type, "] Matrix not found: ", rds_path)
    return(invisible(NULL))
  }

  mats <- readRDS(rds_path)
  mat  <- mats[[mat_type]][mats$accessions, mats$accessions]
  n_na <- sum(is.na(mat)) / 2
  if (n_na > 0)
    message("[", gene, " ", type, " ", mat_type,
            "] Imputing ", n_na, " missing pairs")
  mat <- impute_na(mat)

  message("[", gene, " ", type, " ", mat_type,
          "] Building tree (", nrow(mat), " taxa)...")
  t0   <- proc.time()
  tree <- fastme.bal(as.dist(mat))
  elapsed <- round((proc.time() - t0)[["elapsed"]])
  write.tree(tree, out_nwk)
  message("[", gene, " ", type, " ", mat_type,
          "] Done (", elapsed, "s) → ", basename(out_nwk))
  invisible(NULL)
}

jobs <- unlist(recursive = FALSE,
  unlist(recursive = FALSE, lapply(genes, function(g) {
    lapply(types, function(tp) {
      lapply(c("dN", "dS"), function(mt) {
        list(
          gene     = g,
          type     = tp,
          mat_type = mt,
          out_nwk  = file.path(
            tree_dir,
            paste0(g, "_", tp, "_", mt, "_tree.nwk")
          )
        )
      })
    })
  })))

message("Running ", length(jobs), " tree builds sequentially...")
lapply(jobs, build_one)

message("\nDone.")
