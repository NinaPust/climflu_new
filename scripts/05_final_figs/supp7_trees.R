# scripts/05_final_figs/supp7_trees.R
#
# Supplementary Figure 7: H1 and H3 timetrees (circular) with root-to-tip
# regression plots demonstrating temporal signal.
#
# Layout: 2 columns (H1N1, H3N2), 2 rows (circular tree, root-to-tip scatter).
#
# INPUT:
#   data/trees/timed/{gene}_timetree.nwk
#   data/trees/timed/{gene}_treetime_out/dates.tsv
#   data/trees/timed/{gene}_treetime_out/molecular_clock.txt
#
# OUTPUT: paper/figures/supp7_trees.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ape)
  library(ggtree)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS  <- 9
BS2 <- 10

GENES <- c("H1", "H3")


parse_clock <- function(gene) {
  path <- file.path(PATHS$trees_timed,
                    paste0(gene, "_treetime_out"), "molecular_clock.txt")
  lines <- readLines(path)
  rate  <- as.numeric(sub(".*--rate:\\s*", "", lines[grep("rate", lines)[1]]))
  r2    <- as.numeric(sub(".*--r\\^2:\\s*", "", lines[grep("r\\^2", lines)[1]]))
  list(rate = rate, r2 = r2)
}

make_tree_panel <- function(gene) {
  tree <- read.tree(
    file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
  )

  ggtree(tree, layout = "circular", colour = "black",
         linewidth = 0.1) +
    scale_x_continuous(expand = expansion(mult = 0.002)) +
    scale_y_continuous(expand = expansion(mult = 0.002)) +
    theme(plot.margin = margin(0, 0, 0, 0, "mm"))
}

make_rtt_panel <- function(gene, show_x = TRUE, title = NULL) {
  dates_path <- file.path(PATHS$trees_timed,
                           paste0(gene, "_treetime_out"), "dates.tsv")
  dates <- read_tsv(dates_path, comment = "#",
                    col_names = c("node", "date", "num_date"),
                    show_col_types = FALSE) |>
    filter(!grepl("^NODE_", node))

  clock <- parse_clock(gene)

  div_path <- file.path(PATHS$trees_timed,
                        paste0(gene, "_treetime_out"), "divergence_tree.nexus")
  div_tree <- ape::read.nexus(div_path)
  rtt <- data.frame(
    node     = div_tree$tip.label,
    distance = ape::node.depth.edgelength(div_tree)[seq_len(ape::Ntip(div_tree))]
  )

  plot_df <- dates |>
    inner_join(rtt, by = "node") |>
    filter(is.finite(distance), is.finite(num_date))

  label <- sprintf("rate = %.2e subs/site/yr\nR² = %.2f",
                   clock$rate, clock$r2)

  ggplot(plot_df, aes(x = num_date, y = distance)) +
    geom_point(size = 0.5, alpha = 0.5, colour = viridis(1, begin = 0.5)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5,
                colour = "black", linetype = "dashed") +
    annotate("text", x = Inf, y = -Inf, label = label,
             hjust = 1.05, vjust = -0.3,
             size = (BS - 2) / .pt, colour = "grey30") +
    labs(
      x     = if (show_x) "Collection year" else NULL,
      y     = "Root-to-tip distance (subs/site)",
      title = title
    ) +
    theme_bw(base_size = BS) +
    theme(
      panel.grid.minor     = element_blank(),
      axis.title           = element_text(size = BS2),
      axis.text            = element_text(size = BS - 1),
      axis.text.x          = if (show_x) element_text() else element_blank(),
      axis.ticks.x         = if (show_x) element_line() else element_blank(),
      plot.margin          = margin(2, 4, 2, 4, "mm"),
      legend.position      = "none",
      plot.title           = element_text(size = BS, face = "bold",
                                          hjust = 0, colour = "grey20"),
      plot.title.position  = "plot"
    )
}

p_tree_H1 <- make_tree_panel("H1")
p_tree_H3 <- make_tree_panel("H3")
p_rtt_H1  <- make_rtt_panel("H1", show_x = FALSE, title = "H1N1")
p_rtt_H3  <- make_rtt_panel("H3", show_x = TRUE,  title = "H3N2")

# Row 1: H1 scatter | H1 tree
# Row 2: H3 scatter | H3 tree  (scatter shares x axis with row 1)
fig <- (p_rtt_H1 | p_tree_H1) / (p_rtt_H3 | p_tree_H3) +
  plot_layout(widths = c(1, 5), heights = c(1, 1))

ggsave(file.path(out_dir, "supp7_trees.png"), fig,
       width = 200, height = 240, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp7_trees.pdf"), fig,
       width = 200, height = 240, units = "mm")

message("Written: supp7_trees.png / .pdf")