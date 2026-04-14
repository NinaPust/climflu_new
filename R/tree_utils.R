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
#' @param tree    phylo object.
#' @param min_bl  Minimum allowed terminal branch length (tips shorter removed).
#' @param max_bl  Maximum allowed terminal branch length (tips longer removed).
#' @param round   Round number (integer, for informative messages only).
#' @return Pruned phylo object.
prune_outliers <- function(tree, min_bl = 0, max_bl = Inf, round = 1) {
  # Terminal (tip) branch lengths
  tip_edges   <- which(tree$edge[, 2] <= ape::Ntip(tree))
  tip_lengths <- tree$edge.length[tip_edges]
  tip_labels  <- tree$tip.label[tree$edge[tip_edges, 2]]

  to_drop <- tip_labels[tip_lengths < min_bl | tip_lengths > max_bl]

  if (length(to_drop) == 0) {
    message("  [Round ", round, "] No outlier tips removed (range: ",
            round(min(tip_lengths), 5), " – ", round(max(tip_lengths), 5), ")")
    return(tree)
  }

  message("  [Round ", round, "] Dropping ", length(to_drop), " tips ",
          "(bl < ", min_bl, " or > ", max_bl, "): ",
          paste(head(to_drop, 5), collapse = ", "),
          if (length(to_drop) > 5) "..." else "")

  ape::drop.tip(tree, to_drop)
}
