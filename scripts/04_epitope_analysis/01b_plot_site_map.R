# scripts/04_epitope_analysis/01b_plot_site_map.R
#
# PURPOSE: Linear site map of H1 and H3 HA showing antigenic site
#          positions, HA1/HA2 domains, and signal peptide.
#          Analogous to H5 site map in h5_hostswitch project.
#
# INPUT:   data/epitope/{H1,H3}_site_map.csv
# OUTPUT:  figures/epitope/{H1,H3}_site_map.pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(viridisLite)
})

source(here("config.R"))

out_dir <- here("figures", "epitope")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 5 viridis colours for antigenic sites; ridge hardcoded to viridis yellow
epi_cols <- viridis(5, end = 0.75, option = "viridis")

region_colours <- c(
  signal      = "grey88",
  HA1         = "#C6DCED",
  HA2         = "#6AAFD4",
  epitope_A   = epi_cols[1],
  epitope_B   = epi_cols[2],
  epitope_C   = epi_cols[3],
  epitope_D   = epi_cols[4],
  epitope_E   = epi_cols[5],
  epitope_Sa  = epi_cols[1],
  epitope_Sb  = epi_cols[2],
  epitope_Ca1 = epi_cols[3],
  epitope_Ca2 = epi_cols[4],
  epitope_Cb  = epi_cols[5],
  ridge       = "#FDE725"
)

for (gene in c("H1", "H3")) {
  site_map_path <- file.path(here("data", "epitope"),
                             paste0(gene, "_site_map.csv"))
  if (!file.exists(site_map_path)) {
    warning("Site map not found, skipping ", gene, ": ", site_map_path)
    next
  }

  site_map <- read_csv(site_map_path, show_col_types = FALSE)
  n_aln_cols <- max(site_map$alignment_col)

  # Fill alignment columns not covered by the reference (gap positions in
  # the reference sequence) by propagating the nearest annotated region.
  site_map_full <- tibble(alignment_col = seq_len(n_aln_cols)) |>
    left_join(site_map |> select(alignment_col, region, residue_pos),
              by = "alignment_col") |>
    fill(region, .direction = "downup") |>
    mutate(
      region     = if_else(grepl("_ridge$", region), "ridge", region),
      region     = factor(region, levels = names(region_colours)),
      codon_pos  = ceiling(alignment_col / 3),
      is_epitope = grepl("^epitope_", region) | region == "ridge"
    )

  n_codons <- max(site_map_full$codon_pos)

  blocks <- site_map_full |>
    arrange(codon_pos) |>
    mutate(block = cumsum(c(1, diff(codon_pos) != 1 |
                              region[-1] != region[-n()]))) |>
    group_by(block, region, is_epitope) |>
    summarise(xmin = min(codon_pos) - 0.5,
              xmax = max(codon_pos) + 0.5,
              .groups = "drop")

  p <- ggplot(blocks) +
    geom_rect(aes(xmin = xmin, xmax = xmax,
                  ymin = 0, ymax = 1, fill = region),
              colour = NA) +
    scale_fill_manual(
      values = region_colours,
      breaks = names(region_colours),
      labels = gsub("^epitope_", "", names(region_colours)),
      name   = NULL,
      drop   = TRUE,
      guide  = guide_legend(
        override.aes = list(colour = NA),
        nrow = 2
      )
    ) +
    scale_x_continuous(
      limits = c(0, n_codons),
      expand = c(0, 0),
      breaks = seq(0, n_codons, by = 50)
    ) +
    labs(
      title    = paste0(gene, " HA — antigenic site map"),
      subtitle = paste0("Coloured by domain; x-axis = codon position"),
      x = "Codon position",
      y = NULL
    ) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.y        = element_blank(),
      axis.ticks.y       = element_blank(),
      panel.grid.y       = element_blank(),
      panel.grid.major.x = element_line(colour = "grey90",
                                        linewidth = 0.3),
      panel.grid.minor.x = element_blank(),
      legend.position    = "bottom"
    )

  out_path <- file.path(out_dir, paste0(gene, "_site_map.pdf"))
  ggsave(out_path, p, width = 14, height = 3, dpi = 200)
  message("Saved: ", out_path)
}

message("\nDone.")
