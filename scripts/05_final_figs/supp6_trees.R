# scripts/05_final_figs/supp6_trees.R
#
# Supplementary Figure 7: H1 and H3 timetrees (circular layout), side by side.
#
# INPUT:
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT: paper/figures/supp6_trees.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(ggplot2)
  library(patchwork)
  library(ape)
  library(ggtree)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GENE_TITLES <- c(H1 = "H1", H3 = "H3")

make_tree_panel <- function(gene) {
  tree <- read.tree(
    file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
  )
  ggtree(tree, layout = "circular", colour = "black", linewidth = 0.1) +
    scale_x_continuous(expand = expansion(mult = 0.002)) +
    scale_y_continuous(expand = expansion(mult = 0.002)) +
    labs(title = GENE_TITLES[gene]) +
    theme(
      plot.margin = margin(0, 0, 0, 0, "mm"),
      plot.title  = element_text(size = 9, face = "bold",
                                 hjust = 0, colour = "grey20")
    )
}

p_tree_H1 <- make_tree_panel("H1")
p_tree_H3 <- make_tree_panel("H3")

fig <- p_tree_H1 | p_tree_H3

ggsave(file.path(out_dir, "supp6_trees.png"), fig,
       width = 200, height = 110, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp6_trees.pdf"), fig,
       width = 200, height = 110, units = "mm")

message("Written: supp6_trees.png / .pdf")