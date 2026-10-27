# scripts/05_final_figs/supp8_data_overview.R
#
# Supplementary Figure 4: Dataset overview — distributions of climate and
# demographic variables, and sampling by year and month.
#
#   Row 1: Temperature | Relative humidity | Population density  (density plots)
#   Row 2: Sequences per year | Sequences per month              (bar plots)
#
# INPUT:  data/rates/{gene}_rates.csv
# OUTPUT: paper/figures/supp8_data_overview.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

source(here("revisions", "config.R"))

# Revised figures stay under revisions/paper/ so the submitted ones in
# paper/figures/ are never overwritten.
out_dir <- here("revisions", "paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS  <- 9
BS2 <- 10

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

rates <- lapply(ALL_GENES, function(g) {
  p <- file.path(PATHS$rates, paste0(g, "_rates.csv"))
  if (!file.exists(p)) return(NULL)
  read_csv(p, show_col_types = FALSE) |> mutate(gene = g)
}) |>
  bind_rows() |>
  mutate(
    subtype = case_when(
      gene %in% c("H1", "N1", "PA_H1N1") ~ "H1N1",
      gene %in% c("H3", "N2", "PA_H3N2") ~ "H3N2"
    )
  )

SUBTYPE_COLS <- setNames(viridis(2, end = 0.75), c("H1N1", "H3N2"))

common_theme <- theme_bw(base_size = BS) +
  theme(
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = BS2),
    axis.text          = element_text(size = BS - 1),
    legend.text        = element_text(size = BS - 1),
    legend.key.size    = unit(0.6, "lines"),
    legend.key.spacing.y = unit(0.2, "lines"),
    plot.tag           = element_text(face = "bold.italic", size = BS + 2),
    plot.tag.position  = "topleft"
  )

legend_theme <- theme(
  legend.position      = c(0.05, 0.97),
  legend.justification = c(0, 1),
  legend.direction     = "vertical",
  legend.background    = element_rect(fill = "transparent", colour = NA)
)

clim_df <- rates |>
  filter(!is.na(era5_temp_mean), !is.na(humidity_rh)) |>
  distinct(accession, subtype, era5_temp_mean, humidity_rh, pop_density)

# ── Row 1A: Temperature density ────────────────────────────────────────────────
p_temp <- ggplot(clim_df,
                 aes(x = era5_temp_mean, colour = subtype, fill = subtype)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_colour_manual(values = SUBTYPE_COLS, name = NULL) +
  scale_fill_manual(values = SUBTYPE_COLS, name = NULL) +
  labs(x = "Temperature (°C)", y = "Density", tag = "A") +
  common_theme + legend_theme

# ── Row 1B: Humidity density ───────────────────────────────────────────────────
p_hum <- ggplot(clim_df,
                aes(x = humidity_rh, colour = subtype, fill = subtype)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_colour_manual(values = SUBTYPE_COLS, name = NULL) +
  scale_fill_manual(values = SUBTYPE_COLS, name = NULL) +
  labs(x = "Relative humidity (%)", y = NULL, tag = "B") +
  common_theme + theme(legend.position = "none")

# ── Row 1C: Population density (log10 x-axis) ─────────────────────────────────
p_pop <- ggplot(clim_df |> filter(!is.na(pop_density), pop_density > 0),
                aes(x = pop_density, colour = subtype, fill = subtype)) +
  geom_density(alpha = 0.25, linewidth = 0.7) +
  scale_x_log10(labels = scales::label_comma()) +
  scale_colour_manual(values = SUBTYPE_COLS, name = NULL) +
  scale_fill_manual(values = SUBTYPE_COLS, name = NULL) +
  labs(x = "Population density (per km², log scale)", y = NULL, tag = "C") +
  common_theme + theme(legend.position = "none")

# ── Row 2A: Sequences per year ────────────────────────────────────────────────
gene_labels <- c(
  H1 = "HA (H1N1)", H3 = "HA (H3N2)",
  N1 = "NA (H1N1)", N2 = "NA (H3N2)",
  PA_H1N1 = "PA (H1N1)", PA_H3N2 = "PA (H3N2)"
)
GENE_COLS <- setNames(viridis(6, end = 0.92), gene_labels)

year_df <- rates |>
  filter(!is.na(year)) |>
  distinct(accession, gene, year) |>
  count(gene, year) |>
  mutate(label = factor(gene_labels[gene], levels = gene_labels))

p_year <- ggplot(year_df, aes(x = year, y = n, fill = label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = GENE_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(2009, 2023, 2)) +
  labs(x = "Year", y = "Sequences", tag = "D") +
  common_theme +
  theme(
    panel.grid.major.x   = element_blank(),
    axis.text.x          = element_text(angle = 0, hjust = 0.5),
    legend.position      = c(0.02, 0.97),
    legend.justification = c(0, 1),
    legend.direction     = "vertical",
    legend.background    = element_rect(fill = "transparent", colour = NA)
  )

# ── Row 2B: Sequences per month ───────────────────────────────────────────────
month_df <- rates |>
  filter(!is.na(month)) |>
  distinct(accession, gene, month) |>
  count(gene, month) |>
  mutate(label = factor(gene_labels[gene], levels = gene_labels))

p_month <- ggplot(month_df, aes(x = month, y = n, fill = label)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  scale_fill_manual(values = GENE_COLS, name = NULL) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(x = "Month", y = NULL, tag = "E") +
  common_theme +
  theme(
    panel.grid.major.x = element_blank(),
    axis.text.x        = element_text(angle = 0, hjust = 0.5),
    legend.position    = "none"
  )

# ── Assemble ──────────────────────────────────────────────────────────────────
top_row <- (p_temp | p_hum | p_pop) + plot_layout(nrow = 1)
bot_row <- (p_year | p_month)       + plot_layout(nrow = 1)

fig <- top_row / bot_row +
  plot_layout(heights = c(1, 1))

ggsave(file.path(out_dir, "supp8_data_overview.png"), fig,
       width = 210, height = 160, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp8_data_overview.pdf"), fig,
       width = 210, height = 160, units = "mm")

message("Written: supp8_data_overview.png / .pdf")