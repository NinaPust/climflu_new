# scripts/05_final_figs/supp1_site_map.R
#
# Supplementary Figure 2: Linear site map of H1 and H3 HA showing antigenic
# site positions, HA1/HA2 domains, and signal peptide. H1 on top, H3 below,
# shared legend at the bottom.
#
# INPUT:  data/epitope/{H1,H3}_site_map.csv
# OUTPUT: paper/figures/supp1_site_map.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS <- 9

region_colours <- c(
  signal  = "grey88",
  HA1     = "#C6DCED",
  HA2     = "#6AAFD4",
  epitope = "#7B2D8B",
  ridge   = "#FDE725"
)

make_site_map <- function(gene, show_x = FALSE) {
  site_map_path <- file.path(here("data", "epitope"),
                             paste0(gene, "_site_map.csv"))
  site_map <- read_csv(site_map_path, show_col_types = FALSE)
  n_aln_cols <- max(site_map$alignment_col)

  site_map_full <- tibble(alignment_col = seq_len(n_aln_cols)) |>
    left_join(site_map |> select(alignment_col, region, residue_pos),
              by = "alignment_col") |>
    fill(region, .direction = "downup") |>
    mutate(
      region    = if_else(grepl("_ridge$", region), "ridge", region),
      region    = if_else(grepl("^epitope_", region), "epitope", region),
      region    = factor(region, levels = names(region_colours)),
      codon_pos = ceiling(alignment_col / 3)
    )

  n_codons <- max(site_map_full$codon_pos)

  blocks <- site_map_full |>
    arrange(codon_pos) |>
    mutate(block = cumsum(c(1, diff(codon_pos) != 1 |
                              region[-1] != region[-n()]))) |>
    group_by(block, region) |>
    summarise(xmin = min(codon_pos) - 0.5,
              xmax = max(codon_pos) + 0.5,
              .groups = "drop")

  ggplot(blocks) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = 0, ymax = 1, fill = region),
              colour = NA) +
    scale_fill_manual(
      values = region_colours,
      name   = NULL,
      drop   = TRUE,
      guide  = guide_legend(
        override.aes = list(colour = NA),
        nrow = 1
      )
    ) +
    scale_x_continuous(
      limits = c(0, n_codons),
      expand = c(0, 0),
      breaks = seq(0, n_codons, by = 50)
    ) +
    annotate("text", x = n_codons * 0.01, y = 0.5,
             label = gene, hjust = 0, vjust = 0.5,
             size = (BS + 1) / .pt, fontface = "bold", colour = "grey20") +
    labs(x = if (show_x) "Codon position" else NULL, y = NULL) +
    theme_minimal(base_size = BS) +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.x        = if (show_x) element_text() else element_blank(),
      axis.ticks.x       = if (show_x) element_line() else element_blank(),
      panel.grid.y       = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90", linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      legend.position    = "bottom",
      legend.text        = element_text(size = BS - 1),
      plot.margin        = margin(2, 4, 1, 4, "mm")
    )
}

p_H1 <- make_site_map("H1", show_x = FALSE)
p_H3 <- make_site_map("H3", show_x = TRUE)

fig <- (p_H1 / p_H3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "supp1_site_map.png"), fig,
       width = 180, height = 80, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp1_site_map.pdf"), fig,
       width = 180, height = 80, units = "mm")

message("Written: supp1_site_map.png / .pdf")