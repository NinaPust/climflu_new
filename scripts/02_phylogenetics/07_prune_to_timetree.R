# scripts/02_phylogenetics/07_prune_to_timetree.R
#
# PURPOSE: Prune the already-pruned MAPLE trees to retain only tips that
#          survived TreeTime dating (i.e. are present in the timetree).
#          TreeTime silently drops clock outliers and tips with ambiguous
#          dates; this script makes the pruned-tree tip set match exactly.
#
# INPUT:
#   data/trees/pruned/{gene}_pruned.nwk
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   figures/trees/time_pruned_trees.pdf  — A4 multi-page, one gene per page
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(readr)
})

source(here("config.R"))

args      <- commandArgs(trailingOnly = TRUE)
ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
genes     <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}

out_path <- here("figures", "trees", "time_pruned_trees.pdf")
dir.create(here("figures", "trees"), recursive = TRUE, showWarnings = FALSE)

# A4 in inches
A4_W <- 8.27
A4_H <- 11.69

pdf(out_path, width = A4_W, height = A4_H)

for (gene in genes) {
  message("[", gene, "] pruning to timetree tip set...")

  pruned_path <- file.path(PATHS$trees_pruned, paste0(gene, "_pruned.nwk"))
  timed_path  <- file.path(PATHS$trees_timed,  paste0(gene, "_timetree.nwk"))

  if (!file.exists(pruned_path)) { warning("Missing: ", pruned_path); next }
  if (!file.exists(timed_path))  { warning("Missing: ", timed_path);  next }

  pruned <- read.tree(pruned_path)
  timed  <- read.tree(timed_path)

  # Tips in timetree are the survivors of TreeTime clock filtering
  timed_tips  <- timed$tip.label
  pruned_tips <- pruned$tip.label

  n_before  <- length(pruned_tips)
  keep      <- intersect(pruned_tips, timed_tips)
  dropped   <- setdiff(pruned_tips, timed_tips)
  n_after   <- length(keep)

  message("  ", n_before, " -> ", n_after, " tips  (dropped ", length(dropped), ")")

  tree_final <- keep.tip(pruned, keep)

  # Annotate dropped tips (clock outliers) for colour
  tip_idx  <- tree_final$edge[, 2] <= Ntip(tree_final)
  tip_lens <- tree_final$edge.length[tip_idx]

  long_thr <- quantile(tip_lens, 0.999)   # top 0.1% per gene

  tip_df <- data.frame(
    label    = tree_final$tip.label,
    tip_len  = tip_lens[match(tree_final$tip.label,
                              tree_final$tip.label[tree_final$edge[tip_idx, 2]])],
    stringsAsFactors = FALSE
  ) |>
    mutate(category = if_else(tip_len >= long_thr, "long", "normal"))

  p <- ggtree(tree_final, layout = "rectangular",
              linewidth = 0.04, color = "grey70") %<+% tip_df +
    geom_tippoint(aes(color = category), size = 0.25, shape = 16) +
    scale_color_manual(
      values = c(normal = "grey70", long = "#E69F00"),
      labels = c(normal = "Normal", long = paste0("Long tip (>=", signif(long_thr, 2), ")")),
      name   = NULL
    ) +
    labs(
      title    = paste0(gene, " — pruned to timetree tip set"),
      subtitle = paste0(
        n_before, " tips before  |  ", n_after, " after  |  ",
        length(dropped), " dropped by TreeTime"
      ),
      caption = "Branch lengths: MAPLE substitutions/site. Only tips present in timetree are shown."
    ) +
    theme_tree2() +
    theme(
      legend.position  = "bottom",
      legend.text      = element_text(size = 7),
      plot.title       = element_text(size = 11, face = "bold"),
      plot.subtitle    = element_text(size = 8, color = "grey40"),
      plot.caption     = element_text(size = 6, color = "grey50"),
      plot.background  = element_rect(fill = "white", color = NA),
      panel.background = element_rect(fill = "white", color = NA)
    )

  print(p)
}

dev.off()
message("\nSaved: ", out_path)
