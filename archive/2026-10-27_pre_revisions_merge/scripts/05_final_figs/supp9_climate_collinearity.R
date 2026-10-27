# scripts/05_final_figs/supp9_climate_collinearity.R
#
# Supplementary Figure 9: Temperature vs. relative humidity in the H1N1 and
# H3N2 datasets, addressing potential collinearity between climate predictors.
# Sequences shared across gene segments within a subtype are deduplicated by
# accession so each individual appears once per subtype panel.
#
# INPUT:
#   data/climate/{H1,N1,PA_H1N1,H3,N2,PA_H3N2}_climate.csv
#
# OUTPUT: paper/figures/supp9_climate_collinearity.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS <- 9

# ── Load and pool by subtype, deduplicating by accession ─────────────────────

load_subtype <- function(genes, subtype_label) {
  lapply(genes, function(g) {
    read_csv(here(PATHS$climate, paste0(g, "_climate.csv")),
             show_col_types = FALSE) |>
      select(accession, era5_temp_mean, humidity_rh)
  }) |>
    bind_rows() |>
    distinct(accession, .keep_all = TRUE) |>
    mutate(subtype = subtype_label)
}

dat <- bind_rows(
  load_subtype(c("H1", "N1", "PA_H1N1"), "H1N1"),
  load_subtype(c("H3", "N2", "PA_H3N2"), "H3N2")
)

# ── Per-subtype correlations for annotation ───────────────────────────────────

cors <- dat |>
  group_by(subtype) |>
  summarise(r = cor(era5_temp_mean, humidity_rh, use = "complete.obs"),
            .groups = "drop") |>
  mutate(label = paste0(subtype, "\nr = ", round(r, 3)))

# ── Plot ──────────────────────────────────────────────────────────────────────

p <- ggplot(dat, aes(x = era5_temp_mean, y = humidity_rh)) +
  geom_point(alpha = 0.15, size = 0.6, colour = "#404040") +
  geom_smooth(method = "lm", se = TRUE, colour = "red", linewidth = 0.8) +
  geom_label(data = cors, aes(label = label),
             x = -Inf, y = Inf, hjust = -0.1, vjust = 1.1,
             size = BS * 0.28, fontface = "bold",
             linewidth = 0, fill = "white", inherit.aes = FALSE) +
  facet_wrap(~subtype) +
  labs(
    x = "Monthly mean temperature (°C)",
    y = "Relative humidity (%)"
  ) +
  theme_bw(base_size = BS) +
  theme(
    strip.background = element_blank(),
    strip.text = element_blank()
  )

# ── Save ──────────────────────────────────────────────────────────────────────

ggsave(file.path(out_dir, "supp9_climate_collinearity.png"),
       p, width = 6, height = 3, dpi = 300)
ggsave(file.path(out_dir, "supp9_climate_collinearity.pdf"),
       p, width = 6, height = 3)

message("Saved supp9_climate_collinearity.png/.pdf")
