# scripts/05_final_figs/fig1_data_overview.R
#
# Fig 1: Dataset overview
#   A – World map of sampling locations, coloured by collection temperature
#   B – Temperature and humidity distributions by subtype
#   C – Sequences per year per gene
#
# OUTPUT: paper/figures/fig1_data_overview.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(maps)
  library(viridis)
})

source(here("revisions", "config.R"))

# Revised figures stay under revisions/paper/ so the submitted ones in
# paper/figures/ are never overwritten.
out_dir <- here("revisions", "paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

# ── Load all rates files ───────────────────────────────────────────────────────
rates <- lapply(ALL_GENES, function(g) {
  p <- file.path(PATHS$rates, paste0(g, "_rates.csv"))
  if (!file.exists(p)) return(NULL)
  read_csv(p, show_col_types = FALSE) |>
    mutate(gene = g)
}) |> bind_rows() |>
  mutate(
    subtype  = case_when(
      gene %in% c("H1", "N1", "PA_H1N1") ~ "H1N1",
      gene %in% c("H3", "N2", "PA_H3N2") ~ "H3N2"
    ),
    protein  = case_when(
      gene %in% c("H1", "H3")            ~ "HA",
      gene %in% c("N1", "N2")            ~ "NA",
      TRUE                                ~ "PA"
    )
  )

# One row per unique accession (some genes share the same isolate)
map_df <- rates |>
  filter(!is.na(latitude), !is.na(longitude), !is.na(era5_temp_mean)) |>
  group_by(accession) |>
  slice(1) |>
  ungroup()

# ── Panel A: Sampling map ──────────────────────────────────────────────────────
world <- map_data("world")

pa <- ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey88", colour = "grey70", linewidth = 0.15) +
  geom_point(data = map_df |> slice_sample(n = min(6000, nrow(map_df))),
             aes(x = longitude, y = latitude, colour = era5_temp_mean),
             size = 1.2, alpha = 0.5, shape = 16) +
  scale_colour_viridis_c(name = "Temperature (°C)", option = "plasma",
                         guide = guide_colourbar(barwidth = 0.5, barheight = 5,
                                                 title.position = "top",
                                                 title.hjust = 0.5)) +
  coord_fixed(1.3, xlim = c(-170, 180), ylim = c(-55, 80)) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    panel.grid            = element_blank(),
    axis.text             = element_blank(),
    axis.ticks            = element_blank(),
    # vertical legend, left bottom inside map
    legend.position       = c(0.10, 0.25),
    legend.direction      = "vertical",
    legend.background     = element_rect(fill = "transparent", colour = NA),
    legend.title          = element_text(size = 10),
    legend.text           = element_text(size = 9),
    plot.margin           = margin(2, 2, 2, 2)
  )

# ── Panel B: Climate distributions by subtype ──────────────────────────────────
clim_df <- rates |>
  filter(!is.na(era5_temp_mean), !is.na(humidity_rh)) |>
  distinct(accession, subtype, era5_temp_mean, humidity_rh)

pb_temp <- ggplot(clim_df, aes(x = era5_temp_mean, colour = subtype, fill = subtype)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_colour_viridis_d(end = 0.75, name = NULL) +
  scale_fill_viridis_d(end = 0.75, name = NULL) +
  labs(x = "Temperature (°C)", y = "Density") +
  theme_bw(base_size = 13) +
  theme(
    legend.position    = c(0.05, 0.97),
    legend.justification = c(0, 1),
    legend.direction   = "vertical",
    legend.background  = element_rect(fill = "transparent", colour = NA),
    legend.key.size    = unit(0.7, "lines"),
    legend.key.spacing.y   = unit(0.2, "lines"),
    panel.grid.minor   = element_blank(),
    aspect.ratio       = 1
  )

pb_hum <- ggplot(clim_df,
                 aes(x = humidity_rh, colour = subtype, fill = subtype)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_colour_viridis_d(end = 0.75, name = NULL) +
  scale_fill_viridis_d(end = 0.75, name = NULL) +
  labs(x = "Relative humidity (%)", y = NULL) +
  theme_bw(base_size = 13) +
  theme(legend.position  = "none",
        panel.grid.minor = element_blank(),
        aspect.ratio     = 1)

# ── Panel C: Sequences per year per gene ──────────────────────────────────────
counts_df <- rates |>
  filter(!is.na(year)) |>
  count(gene, year) |>
  mutate(
    protein = case_when(
      gene %in% c("H1", "H3")            ~ "HA",
      gene %in% c("N1", "N2")            ~ "NA",
      TRUE                                ~ "PA"
    ),
    subtype = case_when(
      gene %in% c("H1", "N1", "PA_H1N1") ~ "H1N1",
      TRUE                                ~ "H3N2"
    ),
    label = paste0(protein, "\n(", subtype, ")")
  )

pc <- ggplot(counts_df, aes(x = year, y = n, fill = label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_viridis_d(name = NULL, end = 0.92) +
  scale_x_continuous(breaks = seq(2009, 2023, 2)) +
  labs(x = "Year", y = "Sequences") +
  theme_bw(base_size = 13) +
  theme(
    legend.position      = c(0.02, 0.97),
    legend.justification = c(0, 1),
    legend.direction     = "vertical",
    legend.background    = element_rect(fill = "transparent", colour = NA),
    legend.key.size      = unit(0.7, "lines"),
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    aspect.ratio         = 1,
    axis.text.x          = element_text(angle = 45, hjust = 1)
  )

# ── Panel D: Sequences per month (pooled across years) ────────────────────────
month_df <- rates |>
  filter(!is.na(month)) |>
  count(gene, month) |>
  mutate(
    protein = case_when(
      gene %in% c("H1", "H3")            ~ "HA",
      gene %in% c("N1", "N2")            ~ "NA",
      TRUE                                ~ "PA"
    ),
    subtype = case_when(
      gene %in% c("H1", "N1", "PA_H1N1") ~ "H1N1",
      TRUE                                ~ "H3N2"
    ),
    label = paste0(protein, "\n(", subtype, ")")
  )

pd <- ggplot(month_df, aes(x = month, y = n, fill = label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_viridis_d(name = NULL, end = 0.92,
                       guide = guide_legend(nrow = 1)) +
  scale_x_continuous(breaks = 1:12,
                     labels = month.abb) +
  labs(x = "Month", y = NULL) +
  theme_bw(base_size = 13) +
  theme(
    legend.position  = "none",
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    aspect.ratio     = 1,
    axis.text.x      = element_text(angle = 45, hjust = 1)
  )

# ── Assemble: map full-width on top, four square panels side by side on bottom ─
# Tags: A=map, B=temp density, C=year bar (humidity and month bar unlabelled)
pa      <- pa      + labs(tag = "A")
pb_temp <- pb_temp + labs(tag = "B")
pb_hum  <- pb_hum  + labs(tag = "")
pc      <- pc      + labs(tag = "C")
pd      <- pd      + labs(tag = "")

mid_row <- (pb_temp | pb_hum) +
  plot_layout(nrow = 1)

bot_row <- (pc | pd) + plot_layout(nrow = 1)

fig1 <- pa / mid_row / bot_row +
  plot_layout(heights = c(1.1, 1, 1)) &
  theme(plot.tag          = element_text(face = "bold", size = 12),
        plot.tag.position = "topleft")

ggsave(file.path(out_dir, "fig1_data_overview.png"),
       fig1, width = 8.27, height = 11.69, dpi = 200)
ggsave(file.path(out_dir, "fig1_data_overview.pdf"),
       fig1, width = 8.27, height = 11.69)
message("Saved: paper/figures/fig1_data_overview.png / .pdf")
