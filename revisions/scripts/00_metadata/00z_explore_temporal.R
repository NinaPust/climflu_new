# scripts/00_metadata/00z_explore_temporal.R
#
# Exploratory temporal sampling plots.
#
# OUTPUT:
#   figures/explore_month_year_heatmap.png  — heatmap: month × year, faceted by subtype
#   figures/explore_year_bar.png            — bar plots: counts per year and per month
#   figures/explore_options.png             — 4-panel: dot grid, stacked area, ridge, proportional area

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
  library(ggridges)
  library(tidyr)
})

source(here("revisions", "config.R"))

out_dir <- here("revisions", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Data ──────────────────────────────────────────────────────────────────────

H1N1_GENES <- c("H1", "N1", "PA_H1N1")
H3N2_GENES <- c("H3", "N2", "PA_H3N2")
ALL_GENES  <- c(H1N1_GENES, H3N2_GENES)

load_gene <- function(g) {
  path <- file.path(PATHS$metadata, paste0(g, "_metadata.csv"))
  tryCatch(
    read_csv(path, show_col_types = FALSE) |>
      count(year, month) |>
      mutate(gene = g),
    error = function(e) { warning("Could not load: ", g); NULL }
  )
}

raw <- lapply(ALL_GENES, load_gene) |> bind_rows()

subtype_dat <- bind_rows(
  raw |> filter(gene %in% H1N1_GENES) |>
    group_by(year, month) |> summarise(n = sum(n), .groups = "drop") |>
    mutate(subtype = "H1N1"),
  raw |> filter(gene %in% H3N2_GENES) |>
    group_by(year, month) |> summarise(n = sum(n), .groups = "drop") |>
    mutate(subtype = "H3N2")
)

SUBTYPE_COLS <- setNames(viridis(2, end = 0.85), c("H1N1", "H3N2"))
all_years    <- sort(unique(subtype_dat$year))

common_theme <- theme_bw(base_size = 9) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank())

# ── Plot 1: Heatmap ───────────────────────────────────────────────────────────

p_heat <- ggplot(subtype_dat, aes(x = month, y = year, fill = n)) +
  geom_tile() +
  facet_wrap(~subtype, ncol = 2) +
  coord_fixed() +
  scale_fill_viridis_c(name = "Count", trans = "sqrt") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(breaks = all_years) +
  labs(x = "Month", y = "Year") +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank())

ggsave(file.path(out_dir, "explore_month_year_heatmap.png"), p_heat,
       width = 200, height = 140, units = "mm", dpi = 150)

# ── Plot 2: Bar plots (year + month) ─────────────────────────────────────────

by_year  <- subtype_dat |> group_by(year,  subtype) |> summarise(n = sum(n), .groups = "drop")
by_month <- subtype_dat |> group_by(month, subtype) |> summarise(n = sum(n), .groups = "drop")
ylim_bar <- c(0, max(by_year$n, by_month$n) * 1.05)

p_year <- ggplot(by_year, aes(x = year, y = n, fill = subtype)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = SUBTYPE_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(2000, 2024, 2)) +
  scale_y_continuous(limits = ylim_bar) +
  labs(x = "Year", y = "Sequence count") +
  common_theme +
  theme(legend.position = "top", aspect.ratio = 1)

p_month <- ggplot(by_month, aes(x = month, y = n, fill = subtype)) +
  geom_col(position = "dodge", width = 0.7) +
  scale_fill_manual(values = SUBTYPE_COLS, name = NULL) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(limits = ylim_bar) +
  labs(x = "Month", y = NULL) +
  common_theme +
  theme(legend.position = "none", aspect.ratio = 1)

fig_bar <- p_year | p_month
ggsave(file.path(out_dir, "explore_year_bar.png"), fig_bar,
       width = 200, height = 110, units = "mm", dpi = 150)

# ── Plot 3: 4-panel options ───────────────────────────────────────────────────

by_ym      <- subtype_dat |> group_by(year, month) |> summarise(n = sum(n), .groups = "drop")
by_year_s  <- subtype_dat |> group_by(year,  subtype) |> summarise(n = sum(n), .groups = "drop")

# 3a. Dot grid
p_dot <- ggplot(by_ym, aes(x = month, y = year, size = n)) +
  geom_point(colour = viridis(1, begin = 0.4), alpha = 0.7) +
  scale_size_area(max_size = 6, name = "Count") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_y_continuous(breaks = all_years) +
  labs(title = "Dot grid", x = "Month", y = "Year") +
  theme_bw(base_size = 8) +
  theme(panel.grid = element_blank(), legend.position = "right")

# 3b. Stacked area
p_area <- ggplot(by_year_s, aes(x = year, y = n, fill = subtype)) +
  geom_area(alpha = 0.85) +
  scale_fill_manual(values = SUBTYPE_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(2000, 2024, 4)) +
  labs(title = "Stacked area", x = "Year", y = "Count") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

# 3c. Ridge plot (seasonal density per year)
ridge_dat <- subtype_dat |>
  group_by(year, month) |> summarise(n = sum(n), .groups = "drop") |>
  uncount(n)

p_ridge <- ggplot(ridge_dat, aes(x = month, y = factor(year), fill = after_stat(x))) +
  geom_density_ridges_gradient(scale = 1.5, rel_min_height = 0.01, bandwidth = 1.2) +
  scale_fill_viridis_c(name = "Month") +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(title = "Ridge plot", x = "Month", y = "Year") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

# 3d. Proportional stacked area
p_prop <- ggplot(by_year_s, aes(x = year, y = n, fill = subtype)) +
  geom_area(position = "fill", alpha = 0.85) +
  scale_fill_manual(values = SUBTYPE_COLS, name = NULL) +
  scale_x_continuous(breaks = seq(2000, 2024, 4)) +
  scale_y_continuous(labels = scales::percent) +
  labs(title = "Proportional area", x = "Year", y = "Proportion") +
  theme_bw(base_size = 8) +
  theme(panel.grid.minor = element_blank(), legend.position = "top")

fig_opts <- (p_dot | p_area) / (p_ridge | p_prop)
ggsave(file.path(out_dir, "explore_options.png"), fig_opts,
       width = 220, height = 200, units = "mm", dpi = 150)

message("Written: explore_month_year_heatmap.png, explore_year_bar.png, explore_options.png")
