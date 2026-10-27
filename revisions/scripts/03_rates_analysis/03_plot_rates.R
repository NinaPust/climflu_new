# scripts/03_rates_analysis/03_plot_rates.R
#
# PURPOSE: Exploratory plots of per-tip evolutionary rates and predictors.
#          All panels saved to a single PDF.
#
# INPUT:   data/rates/{gene}_rates.csv
# OUTPUT:  figures/rates_exploration.pdf
#
# USAGE:   Rscript scripts/03_rates_analysis/03_plot_rates.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(tidyr)
  library(scales)
})

source(here("revisions", "config.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

GENE_LABELS <- c(
  H1 = "H1 (HA / H1N1)",  H3 = "H3 (HA / H3N2)",
  N1 = "N1 (NA / H1N1)",  N2 = "N2 (NA / H3N2)",
  PA_H1N1 = "PA (H1N1)",  PA_H3N2 = "PA (H3N2)"
)

GENE_COLOURS <- setNames(
  viridisLite::viridis(length(ALL_GENES), end = 0.92),
  ALL_GENES
)

# ── Load rates data ────────────────────────────────────────────────────────────
dat <- lapply(ALL_GENES, function(g) {
  path <- file.path(PATHS$rates, paste0(g, "_rates.csv"))
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  read_csv(path, show_col_types = FALSE) |>
    mutate(gene   = g,
           label  = GENE_LABELS[g],
           subtype = ifelse(g %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2"),
           protein = case_when(
             g %in% c("H1", "H3")            ~ "HA",
             g %in% c("N1", "N2")            ~ "NA (neuraminidase)",
             g %in% c("PA_H1N1", "PA_H3N2") ~ "PA (control)"
           ))
}) |> bind_rows()

out_dir <- here("revisions", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
pdf_path <- file.path(out_dir, "rates_exploration.pdf")
pdf(pdf_path, width = 11, height = 7)

# ── 1. Log rate distribution per gene ─────────────────────────────────────────
p1 <- dat |>
  filter(is.finite(log_rate)) |>
  ggplot(aes(x = log_rate, fill = gene, colour = gene)) +
  geom_density(alpha = 0.35, linewidth = 0.6) +
  scale_fill_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
  scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
  labs(title = "Distribution of log terminal branch rate by gene",
       x = "log(rate)  [subst / site / yr]", y = "Density") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p1)

# ── 3. Rate over time ──────────────────────────────────────────────────────────
p3 <- dat |>
  filter(is.finite(log_rate)) |>
  ggplot(aes(x = year, y = log_rate, colour = gene)) +
  geom_point(alpha = 0.15, size = 0.7) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 1, span = 0.6) +
  facet_wrap(~subtype) +
  scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
  labs(title = "Terminal branch rate over time",
       x = "Year", y = "log(rate)  [subst / site / yr]") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p3)

# ── 4. Rate vs temperature ─────────────────────────────────────────────────────
p4 <- dat |>
  filter(is.finite(log_rate), is.finite(era5_temp_mean)) |>
  ggplot(aes(x = era5_temp_mean, y = log_rate, colour = gene)) +
  geom_point(alpha = 0.12, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
  facet_wrap(~protein, scales = "free_y") +
  scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
  labs(title = "Rate vs ERA5 temperature",
       x = "Monthly mean temperature (°C)", y = "log(rate)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p4)

# ── 5. Rate vs relative humidity ──────────────────────────────────────────────
p5 <- dat |>
  filter(is.finite(log_rate), is.finite(humidity_rh)) |>
  ggplot(aes(x = humidity_rh, y = log_rate, colour = gene)) +
  geom_point(alpha = 0.12, size = 0.8) +
  geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
  facet_wrap(~protein, scales = "free_y") +
  scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
  labs(title = "Rate vs relative humidity",
       x = "Relative humidity (%)", y = "log(rate)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p5)

# ── 6. Rate by month (seasonal pattern) ───────────────────────────────────────
p6 <- dat |>
  filter(is.finite(log_rate)) |>
  mutate(month = factor(month, levels = 1:12,
                        labels = month.abb)) |>
  ggplot(aes(x = month, y = log_rate, fill = gene)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3,
               linewidth = 0.4, width = 0.7) +
  facet_wrap(~label, scales = "free_y", ncol = 3) +
  scale_fill_manual(values = GENE_COLOURS, guide = "none") +
  labs(title = "Rate by month of sample collection",
       x = "Month", y = "log(rate)") +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
        strip.text  = element_text(size = 9))
print(p6)

# ── 7. Season cosine vs rate ───────────────────────────────────────────────────
if ("season_cos" %in% names(dat)) {
  p7 <- dat |>
    filter(is.finite(log_rate), is.finite(season_cos)) |>
    ggplot(aes(x = season_cos, y = log_rate, colour = gene)) +
    geom_point(alpha = 0.12, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    facet_wrap(~protein, scales = "free_y") +
    scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
    scale_x_continuous(
      breaks = c(-1, 0, 1),
      labels = c("-1\n(off-season)", "0", "+1\n(peak flu season)")
    ) +
    labs(title = "Rate vs hemisphere-adjusted flu season",
         subtitle = "season_cos = cos(2*pi*adj_month/12); SH months shifted +6 to align peak",
         x = "Season cosine", y = "log(rate)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  print(p7)
}

# ── 8. Temperature vs humidity (predictor correlation) ────────────────────────
p8 <- dat |>
  filter(is.finite(era5_temp_mean), is.finite(humidity_rh)) |>
  ggplot(aes(x = era5_temp_mean, y = humidity_rh, colour = gene)) +
  geom_point(alpha = 0.10, size = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
  facet_wrap(~subtype) +
  scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
  labs(title = "Predictor correlation: temperature vs humidity",
       x = "Monthly mean temperature (°C)", y = "Relative humidity (%)") +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")
print(p8)

# ── 9. Geographic distribution of rates ───────────────────────────────────────
world <- map_data("world")
pts9 <- dat |>
  filter(is.finite(log_rate), is.finite(latitude), is.finite(longitude))
p9 <- ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey92", colour = "grey75", linewidth = 0.15) +
  geom_point(data = pts9,
             aes(x = longitude, y = latitude, colour = log_rate),
             alpha = 0.4, size = 0.8) +
  facet_wrap(~label, ncol = 3) +
  scale_colour_viridis_c(name = "log(rate)", option = "plasma") +
  coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-60, 75)) +
  labs(title = "Geographic distribution of terminal branch rates") +
  theme_bw(base_size = 10) +
  theme(axis.text       = element_blank(),
        axis.ticks      = element_blank(),
        axis.title      = element_blank(),
        strip.text      = element_text(size = 8),
        legend.position = "bottom",
        panel.grid      = element_blank())
print(p9)

# ── 10. Sample counts per year x gene ─────────────────────────────────────────
p10 <- dat |>
  count(gene, year) |>
  ggplot(aes(x = year, y = n, fill = gene)) +
  geom_col(width = 0.8) +
  facet_wrap(~gene, labeller = as_labeller(GENE_LABELS), scales = "free_y") +
  scale_fill_manual(values = GENE_COLOURS, guide = "none") +
  labs(title = "Sample counts per year by gene",
       x = "Year", y = "Number of tips") +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(size = 9))
print(p10)

# ── Omega section (only if log_omega data are present) ────────────────────────
has_omega_data <- dat |> filter(is.finite(log_omega)) |> nrow() > 0

if (has_omega_data) {

  # ── O1. Log omega distribution per gene ───────────────────────────────────────
  pO1 <- dat |>
    filter(is.finite(log_omega)) |>
    ggplot(aes(x = log_omega, fill = gene, colour = gene)) +
    geom_density(alpha = 0.35, linewidth = 0.6) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey40") +
    scale_fill_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
    scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
    labs(title = "Distribution of log(dN/dS) by gene",
         subtitle = "Dashed line at log(dN/dS) = 0, i.e. dN/dS = 1 (neutral evolution)",
         x = "log(dN/dS)", y = "Density") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  print(pO1)

  # ── O2. Omega over time ────────────────────────────────────────────────────────
  pO2 <- dat |>
    filter(is.finite(log_omega)) |>
    ggplot(aes(x = year, y = log_omega, colour = gene)) +
    geom_point(alpha = 0.15, size = 0.7) +
    geom_smooth(method = "loess", se = FALSE, linewidth = 1, span = 0.6) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~subtype) +
    scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
    labs(title = "log(dN/dS) over time",
         x = "Year", y = "log(dN/dS)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  print(pO2)

  # ── O3. Omega vs temperature ───────────────────────────────────────────────────
  pO3 <- dat |>
    filter(is.finite(log_omega), is.finite(era5_temp_mean)) |>
    ggplot(aes(x = era5_temp_mean, y = log_omega, colour = gene)) +
    geom_point(alpha = 0.12, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~protein, scales = "free_y") +
    scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
    labs(title = "log(dN/dS) vs ERA5 temperature",
         x = "Monthly mean temperature (°C)", y = "log(dN/dS)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  print(pO3)

  # ── O4. Omega vs relative humidity ────────────────────────────────────────────
  pO4 <- dat |>
    filter(is.finite(log_omega), is.finite(humidity_rh)) |>
    ggplot(aes(x = humidity_rh, y = log_omega, colour = gene)) +
    geom_point(alpha = 0.12, size = 0.8) +
    geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~protein, scales = "free_y") +
    scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
    labs(title = "log(dN/dS) vs relative humidity",
         x = "Relative humidity (%)", y = "log(dN/dS)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")
  print(pO4)

  # ── O5. Omega by month ────────────────────────────────────────────────────────
  pO5 <- dat |>
    filter(is.finite(log_omega)) |>
    mutate(month = factor(month, levels = 1:12, labels = month.abb)) |>
    ggplot(aes(x = month, y = log_omega, fill = gene)) +
    geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3,
                 linewidth = 0.4, width = 0.7) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    facet_wrap(~label, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = GENE_COLOURS, guide = "none") +
    labs(title = "log(dN/dS) by month of sample collection",
         x = "Month", y = "log(dN/dS)") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          strip.text  = element_text(size = 9))
  print(pO5)

  # ── O6. Season cosine vs omega ─────────────────────────────────────────────────
  if ("season_cos" %in% names(dat)) {
    pO6 <- dat |>
      filter(is.finite(log_omega), is.finite(season_cos)) |>
      ggplot(aes(x = season_cos, y = log_omega, colour = gene)) +
      geom_point(alpha = 0.12, size = 0.8) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 0.9, alpha = 0.15) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      facet_wrap(~protein, scales = "free_y") +
      scale_colour_manual(values = GENE_COLOURS, labels = GENE_LABELS, name = NULL) +
      scale_x_continuous(
        breaks = c(-1, 0, 1),
        labels = c("-1\n(off-season)", "0", "+1\n(peak flu season)")
      ) +
      labs(title = "log(dN/dS) vs hemisphere-adjusted flu season",
           x = "Season cosine", y = "log(dN/dS)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")
    print(pO6)
  }

  # ── O7. Geographic distribution of omega ──────────────────────────────────────
  ptsO7 <- dat |>
    filter(is.finite(log_omega), is.finite(latitude), is.finite(longitude))
  pO7 <- ggplot() +
    geom_polygon(data = world,
                 aes(x = long, y = lat, group = group),
                 fill = "grey92", colour = "grey75", linewidth = 0.15) +
    geom_point(data = ptsO7,
               aes(x = longitude, y = latitude, colour = log_omega),
               alpha = 0.4, size = 0.8) +
    facet_wrap(~label, ncol = 3) +
    scale_colour_viridis_c(name = "log(dN/dS)", option = "plasma") +
    coord_fixed(ratio = 1.3, xlim = c(-180, 180), ylim = c(-60, 75)) +
    labs(title = "Geographic distribution of log(dN/dS)") +
    theme_bw(base_size = 10) +
    theme(axis.text       = element_blank(),
          axis.ticks      = element_blank(),
          axis.title      = element_blank(),
          strip.text      = element_text(size = 8),
          legend.position = "bottom",
          panel.grid      = element_blank())
  print(pO7)

} # end has_omega_data

dev.off()
message("Saved: ", pdf_path)
