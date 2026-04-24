# R/tree_utils.R
# Utilities for tree loading and outlier pruning.

#' Load a Newick tree for a given gene segment.
#'
#' @param gene      Gene segment name (e.g. "H1").
#' @param type      Tree type: "tree" | "timetree" | "dN_tree" | "dS_tree"
#' @param trees_dir Path to trees directory (defaults to PATHS$trees from config).
#' @return phylo object.
load_tree <- function(gene, type = "timetree", trees_dir = PATHS$trees) {
  path <- file.path(trees_dir, paste0(gene, "_", type, ".nwk"))
  if (!file.exists(path)) stop("Tree file not found: ", path)
  ape::read.tree(path)
}


#' Prune tips with outlier branch lengths from a phylo tree.
#'
#' Short branches (< min_bl) are removed with a fixed cutoff — these are
#' near-identical duplicate sequences and the threshold is unambiguous.
#'
#' Long branches are detected statistically: the log-transformed tip lengths
#' follow an approximately normal distribution for the bulk of sequences, so
#' robust z-scores (using median and MAD) identify outliers without assuming
#' a fixed cutoff. A tip is flagged if:
#'   (log(bl) - median(log(bl))) / MAD(log(bl)) > k_long
#'
#' @param tree    phylo object.
#' @param min_bl  Fixed lower cutoff: tips shorter than this are removed.
#' @param k_long  MAD multiplier for long-branch detection (default 5).
#'                Higher = fewer tips removed. Tune if too aggressive/lenient.
#' @param round   Round number (integer, for informative messages only).
#' @return Pruned phylo object.
prune_outliers <- function(tree, min_bl = 0, k_long = 5, round = 1) {
  tip_edges   <- which(tree$edge[, 2] <= ape::Ntip(tree))
  tip_lengths <- tree$edge.length[tip_edges]
  tip_labels  <- tree$tip.label[tree$edge[tip_edges, 2]]

  # Fixed lower cutoff — removes identical/near-identical sequences
  too_short <- tip_lengths < min_bl

  # Statistical upper cutoff on log scale — robust to multiple outliers
  # Small offset avoids log(0) for exactly-zero branches
  log_bl  <- log(tip_lengths + 1e-9)
  med_log <- median(log_bl)
  mad_log <- mad(log_bl)
  max_bl  <- exp(med_log + k_long * mad_log)
  too_long <- tip_lengths > max_bl

  to_drop <- tip_labels[too_short | too_long]

  message("  [Round ", round, "] Long-branch threshold (k=", k_long, "): ",
          round(max_bl, 5),
          "  Flagged: ", sum(too_short), " short, ", sum(too_long), " long")

  if (length(to_drop) == 0) {
    message("  [Round ", round, "] No outlier tips removed")
    return(tree)
  }

  message("  [Round ", round, "] Dropping ", length(to_drop), " tips: ",
          paste(head(to_drop, 5), collapse = ", "),
          if (length(to_drop) > 5) "..." else "")

  ape::drop.tip(tree, to_drop)
}
