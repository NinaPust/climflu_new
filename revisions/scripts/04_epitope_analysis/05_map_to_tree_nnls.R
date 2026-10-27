# scripts/04_epitope_analysis/05_map_to_tree_nnls.R
#
# PURPOSE: Fit dN and dS branch lengths onto the timetree topology using
#          non-negative least squares (nnls.tree / phangorn), producing
#          trees with the same topology as the timetree but dN/dS-scaled
#          branch lengths. All 8 jobs run in parallel via mclapply.
#          Use 05b_map_to_tree_fastme.R as a faster alternative; whichever
#          finishes first feeds into 06_compute_epitope_rates.R.
#
# INPUT:
#   data/epitope/{gene}_{type}_dNdS_matrix.rds
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/epitope/trees/{gene}_{type}_{dN|dS}_tree.nwk
#
# FLAGS:
#   --gene X   run a single gene (default: H1 H3)
#   --type X   epitope|nonepitope (default: both)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(phangorn)
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

fit_one <- function(job) {
  gene     <- job$gene
  type     <- job$type
  mat_type <- job$mat_type
  out_nwk  <- job$out_nwk

  if (file.exists(out_nwk)) {
    message("[", gene, " ", type, " ", mat_type, "] Already exists, skipping.")
    return(invisible(NULL))
  }

  rds_path  <- file.path(epi_dir,
                         paste0(gene, "_", type, "_dNdS_matrix.rds"))
  tree_path <- file.path(PATHS$trees_timed,
                         paste0(gene, "_timetree.nwk"))

  if (!file.exists(rds_path)) {
    warning("[", gene, " ", type, "] Matrix not found: ", rds_path)
    return(invisible(NULL))
  }
  if (!file.exists(tree_path)) {
    warning("[", gene, " ", type, "] Timetree not found: ", tree_path)
    return(invisible(NULL))
  }

  mats <- readRDS(rds_path)
  tree <- read.tree(tree_path)

  common <- intersect(mats$accessions, tree$tip.label)
  tree   <- keep.tip(tree, common)

  # Drop zero-length tips (nnls.tree fails on degenerate branches)
  tip_idx   <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips <- tree$tip.label[tree$edge.length[tip_idx] == 0]
  if (length(zero_tips) > 0) {
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  mat <- mats[[mat_type]][common, common]
  n_na <- sum(is.na(mat)) / 2
  if (n_na > 0)
    message("[", gene, " ", type, " ", mat_type,
            "] Imputing ", n_na, " missing pairs")
  mat <- impute_na(mat)

  message("[", gene, " ", type, " ", mat_type,
          "] Fitting nnls.tree (", length(common), " taxa)...")
  t0 <- proc.time()
  fitted <- tryCatch(
    nnls.tree(as.dist(mat), tree, rooted = FALSE),
    error = function(e) {
      message("[", gene, " ", type, " ", mat_type,
              "] nnls.tree error: ", e$message)
      NULL
    }
  )
  elapsed <- round((proc.time() - t0)[["elapsed"]])

  if (is.null(fitted)) return(invisible(NULL))

  write.tree(fitted, out_nwk)
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

message("Running ", length(jobs), " nnls.tree fits sequentially...")
lapply(jobs, fit_one)

message("\nDone.")
