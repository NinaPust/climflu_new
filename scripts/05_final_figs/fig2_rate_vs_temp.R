# scripts/05_final_figs/fig2_rate_vs_temp.R
#
# Fig 2: Evolutionary rate vs. collection temperature
#   Scatter of log terminal branch rate against ERA5 temperature,
#   faceted by subtype (H1N1 / H3N2) and protein (HA / NA).
#   Points coloured by year; LOESS smoother with 95% CI.
#   PA panel included to show null pattern (negative control).
#
# OUTPUT: paper/figures/fig2_rate_vs_temp.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

# ── Load data ─────────────────────────────────────────────────────────────────
rates <- lapply(ALL_GENES, function(g) {
  p <- file.path(PATHS$rates, paste0(g, "_rates.csv"))
  if (!file.exists(p)) return(NULL)
  read_csv(p, show_col_types = FALSE) |> mutate(gene = g)
}) |> bind_rows() |>
  filter(!is.na(era5_temp_mean), !is.na(log_rate)) |>
  mutate(
    subtype = case_when(
      gene %in% c("H1", "N1", "PA_H1N1") ~ "H1N1",
      TRUE                                ~ "H3N2"
    ),
    protein = factor(
      case_when(
        gene %in% c("H1", "H3")           ~ "HA",
        gene %in% c("N1", "N2")           ~ "NA",
        TRUE                               ~ "PA (control)"
      ),
      levels = c("HA", "NA", "PA (control)")
    ),
    facet_label = paste0(protein, " — ", subtype)
  )

# Ordered facets: HA H1N1, HA H3N2, NA H1N1, NA H3N2, PA H1N1, PA H3N2
rates <- rates |>
  mutate(facet_label = factor(facet_label, levels = c(
    "HA — H1N1", "HA — H3N2",
    "NA — H1N1", "NA — H3N2",
    "PA (control) — H1N1", "PA (control) — H3N2"
  )))

# Subsample for plot clarity (keep all for smoother, thin points)
set.seed(42)
plot_pts <- rates |>
  group_by(gene) |>
  mutate(row_id = row_number(), total = n()) |>
  filter(row_id <= pmin(total, 2000)) |>
  select(-row_id, -total) |>
  ungroup()

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(rates, aes(x = era5_temp_mean, y = log_rate)) +

  # LOESS smoother on full data
  geom_smooth(method = "loess", span = 0.6,
              colour = "firebrick", fill = "firebrick",
              alpha = 0.15, linewidth = 0.9, se = TRUE) +

  # Thinned scatter points
  geom_point(data = plot_pts,
             aes(colour = year), size = 0.5, alpha = 0.35, shape = 16) +

  scale_colour_viridis_c(name = "Year", option = "viridis",
                         guide = guide_colourbar(barwidth = 6, barheight = 0.4,
                                                 title.position = "top")) +

  facet_wrap(~ facet_label, ncol = 2, scales = "free_y") +

  labs(
    x = "Collection temperature (°C)",
    y = "log(evolutionary rate)"
  ) +

  theme_bw(base_size = 11) +
  theme(
    strip.text       = element_text(face = "bold", size = 10),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave(file.path(out_dir, "fig2_rate_vs_temp.png"),
       p, width = 9, height = 11, dpi = 200)
message("Saved: paper/figures/fig2_rate_vs_temp.png")
