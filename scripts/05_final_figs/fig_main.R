# scripts/05_final_figs/fig_main.R
#
# Main composite publication figure — A4 width, 4 rows:
#
#   Row 1: Map (black dots) | Forest plot — all genes, H1N1 + H3N2
#   Row 2: Scatter log(rate) vs temperature, H1 | Multi-level forest, H1
#   Row 3: Jackknife ordered strip — H1
#   Row 4: Jackknife ordered strip — H3
#
# References: fig1, fig2, fig3, fig4, fig_jk_b scripts in this directory.
# OUTPUT: paper/figures/fig_main.png / fig_main.pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(cowplot)
  library(maps)
  library(viridis)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS  <- 10   # base font size
BS2 <- 11   # slightly larger for axis titles

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

common_theme <- theme_bw(base_size = BS) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = BS + 1),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = BS2),
    axis.text          = element_text(size = BS - 1),
    legend.text        = element_text(size = BS - 1),
    legend.key.size    = unit(0.6, "lines"),
    plot.tag           = element_text(face = "bold.italic", size = BS + 2),
    plot.tag.position  = "topleft"
  )

# ══════════════════════════════════════════════════════════════════════════════
# ROW 1-A: Map with black dots
# ══════════════════════════════════════════════════════════════════════════════
rates_all <- lapply(ALL_GENES, function(g) {
  p <- file.path(PATHS$rates, paste0(g, "_rates.csv"))
  if (!file.exists(p)) return(NULL)
  read_csv(p, show_col_types = FALSE) |> mutate(gene = g)
}) |> bind_rows()

map_df <- rates_all |>
  filter(!is.na(latitude), !is.na(longitude)) |>
  group_by(accession) |> slice(1) |> ungroup()

set.seed(42)
world <- map_data("world")

p_map <- ggplot() +
  geom_polygon(data = world,
               aes(x = long, y = lat, group = group),
               fill = "grey88", colour = "grey70", linewidth = 0.12) +
  geom_point(data = map_df |> slice_sample(n = min(5000, nrow(map_df))),
             aes(x = longitude, y = latitude),
             size = 0.5, alpha = 0.35, colour = "black", shape = 16) +
  coord_fixed(1.3, xlim = c(-170, 180), ylim = c(-55, 80)) +
  labs(x = NULL, y = NULL, tag = "A") +
  common_theme +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank(),
        panel.grid = element_blank())

# ══════════════════════════════════════════════════════════════════════════════
# ROW 1-B: Forest plot — log(branch rate) ~ all predictors, all genes
# ══════════════════════════════════════════════════════════════════════════════
KEEP_TERMS <- c("scaled_temp", "scaled_humidity",
                "scaled_sampling_density", "scaled_delta_temp",
                "scaled_delta_humidity")
TERM_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative\nhumidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative\nhumidity",
  scaled_sampling_density = "Sampling\ndensity"
)
PRED_COLOURS <- setNames(viridis(length(TERM_LABELS), end = 0.92),
                         TERM_LABELS)

res_main <- read_csv(
  file.path(PATHS$rates_phylolm, "all_phylolm_results.csv"),
  show_col_types = FALSE
) |>
  filter(response == "log_rate", term %in% KEEP_TERMS) |>
  mutate(
    ci_lo     = estimate - 1.96 * std_error,
    ci_hi     = estimate + 1.96 * std_error,
    sig_fdr   = !is.na(p.fdr) & p.fdr < 0.05,
    protein   = factor(
      case_when(
        gene %in% c("H1", "H3")           ~ "HA",
        gene %in% c("N1", "N2")           ~ "NA",
        TRUE                               ~ "PA"
      ), levels = c("HA", "NA", "PA")
    ),
    subtype   = ifelse(gene %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2"),
    predictor = factor(TERM_LABELS[term], levels = TERM_LABELS),
    pt_fill   = ifelse(sig_fdr,
                       PRED_COLOURS[as.character(predictor)],
                       "white")
  )

dodge <- position_dodge(width = 0.5)

forest_labels <- data.frame(subtype = c("H1N1", "H3N2"))

p_forest <- ggplot(res_main,
                   aes(x = estimate, y = protein,
                       colour = predictor, group = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = pt_fill),
             shape = 21, size = 2.2, stroke = 0.7, position = dodge) +
  geom_text(data = forest_labels,
            aes(x = -Inf, y = Inf, label = subtype),
            hjust = -0.1, vjust = 1.4, inherit.aes = FALSE,
            size = (BS - 1) / .pt, colour = "grey20", fontface = "bold") +
  scale_colour_manual(values = PRED_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev) +
  facet_wrap(~ subtype, ncol = 1) +
  labs(x = "Regression coefficient", y = NULL, tag = "C") +
  common_theme +
  theme(legend.position      = c(-0.01, 0.01),
        legend.justification = c(0, 0),
        legend.direction     = "vertical",
        legend.background    = element_rect(fill = "transparent", colour = NA),
        legend.text          = element_text(size = BS - 2),
        legend.key.spacing.y = unit(0.2, "lines"),
        panel.grid.major.x   = element_blank(),
        aspect.ratio         = 0.5,
        strip.text           = element_blank(),
        strip.background     = element_blank()) +
  guides(colour = guide_legend(nrow = 5))

# ══════════════════════════════════════════════════════════════════════════════
# ROW 2-A: Scatter — log(rate) vs temperature, H1 only
# ══════════════════════════════════════════════════════════════════════════════
set.seed(42)
pts_scatter <- rates_all |>
  filter(gene %in% c("H1", "H3"), !is.na(era5_temp_mean), !is.na(log_rate)) |>
  mutate(subtype = ifelse(gene == "H1", "H1N1", "H3N2")) |>
  group_by(subtype) |>
  mutate(row_id = row_number(), total = n()) |>
  filter(row_id <= pmin(total, 2000)) |>
  select(-row_id, -total) |>
  ungroup()

SCATTER_COL <- setNames(viridis(2, end = 0.75), c("H1N1", "H3N2"))

make_scatter <- function(dat, tag_letter, show_y = TRUE, show_x = TRUE) {
  lbl <- unique(dat$subtype)
  p <- ggplot(dat, aes(x = era5_temp_mean, y = log_rate)) +
    geom_point(size = 0.5, alpha = 0.3, shape = 16,
               colour = SCATTER_COL[lbl]) +
    annotate("text", x = -Inf, y = Inf, label = lbl,
             hjust = -0.1, vjust = 1.4,
             size = (BS - 1) / .pt, colour = "grey20", fontface = "bold") +
    labs(x = if (show_x) "Temperature (°C)" else NULL,
         y = if (show_y) "log(substitution rate)" else NULL,
         tag = tag_letter) +
    common_theme +
    theme(aspect.ratio = 0.8)
  if (!show_y) p <- p + theme(axis.text.y = element_blank(),
                               axis.ticks.y = element_blank())
  p
}

p_scatter_h1 <- make_scatter(filter(pts_scatter, subtype == "H1N1"), "B",
                              show_y = TRUE,  show_x = TRUE)
p_scatter_h3 <- make_scatter(filter(pts_scatter, subtype == "H3N2"), "",
                              show_y = FALSE, show_x = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# ROW 2-B: Multi-level forest — temperature coefficient, H1 only
# ══════════════════════════════════════════════════════════════════════════════
ml_main <- read_csv(
  file.path(PATHS$rates_phylolm, "all_phylolm_results.csv"),
  show_col_types = FALSE
) |>
  filter(response == "log_rate", term == "scaled_temp", gene %in% c("H1", "H3")) |>
  mutate(region = "Terminal\nbranch rate")

ml_ns <- read_csv(
  file.path(PATHS$rates_phylolm, "all_phylolm_ns_results.csv"),
  show_col_types = FALSE
) |>
  filter(term == "scaled_temp", gene %in% c("H1", "H3")) |>
  mutate(region = case_when(
    response == "log_full_ns_rate"    ~ "Whole-gene\nNS rate",
    response == "log_full_s_rate"     ~ "Whole-gene\nS rate",
    response == "log_epitope_ns_rate" ~ "Epitope\nNS rate",
    response == "log_ridge_ns_rate"   ~ "Ridge\nNS rate",
    response == "log_HA2_ns_rate"     ~ "HA2\nNS rate",
    TRUE                              ~ NA_character_
  )) |>
  filter(!is.na(region))

ml_combined <- bind_rows(ml_main, ml_ns) |>
  mutate(
    subtype   = ifelse(gene == "H1", "H1N1", "H3N2"),
    region = factor(region, levels = c(
      "Terminal\nbranch rate", "Whole-gene\nNS rate",
      "Whole-gene\nS rate", "HA2\nNS rate", "Epitope\nNS rate",
      "Ridge\nNS rate"
    )),
    ci_lo     = estimate - 1.96 * std_error,
    ci_hi     = estimate + 1.96 * std_error,
    sig_fdr   = !is.na(p.fdr) & p.fdr < 0.05,
    rate_type = factor(case_when(
      region == "Terminal\nbranch rate"  ~ "Substitution\nrate",
      grepl("NS", region)                ~ "Non-synonymous\nrate",
      TRUE                               ~ "Synonymous\nrate"
    ), levels = c("Substitution\nrate", "Non-synonymous\nrate", "Synonymous\nrate"))
  )

TYPE_COLOURS <- setNames(viridis(3, end = 0.85),
                         c("Substitution\nrate", "Non-synonymous\nrate", "Synonymous\nrate"))

ml_labels <- data.frame(subtype = c("H1N1", "H3N2"))

p_multilevel <- ggplot(ml_combined,
                       aes(x = estimate, y = region, colour = rate_type)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.9) +
  geom_point(aes(fill = ifelse(sig_fdr,
                               TYPE_COLOURS[as.character(rate_type)],
                               "white")),
             shape = 21, size = 3, stroke = 0.9) +
  geom_text(data = ml_labels,
            aes(x = -Inf, y = Inf, label = subtype),
            hjust = -0.1, vjust = 1.4, inherit.aes = FALSE,
            size = (BS - 1) / .pt, colour = "grey20", fontface = "bold") +
  scale_colour_manual(values = TYPE_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev, expand = expansion(add = 0.4)) +
  facet_wrap(~ subtype, ncol = 1) +
  geom_text(data = data.frame(subtype = "H3N2"),
            aes(x = Inf, y = -Inf,
                label = "● FDR < 0.05   ○ n.s."),
            hjust = 1.05, vjust = -0.4, inherit.aes = FALSE,
            size = (BS - 2) / .pt, colour = "grey30") +
  labs(x = "Temperature coefficient", y = NULL, tag = "D") +
  common_theme +
  theme(legend.position    = "none",
        panel.grid.major.x = element_blank(),
        aspect.ratio       = 0.5,
        strip.text         = element_blank(),
        strip.background   = element_blank())

# ══════════════════════════════════════════════════════════════════════════════
# ROWS 3 & 4: Jackknife ordered strip — H1 and H3 separately
# ══════════════════════════════════════════════════════════════════════════════
SIG_LEVELS  <- c("ns", "p < 0.05", "p < 0.01", "p < 0.001")
SIG_COLOURS <- setNames(viridis(4, end = 0.92), SIG_LEVELS)

jk_all <- bind_rows(
  read_csv(file.path(PATHS$rates, "jackknife", "H1",
                     "jackknife_results.csv"), show_col_types = FALSE),
  read_csv(file.path(PATHS$rates, "jackknife", "H3",
                     "jackknife_results.csv"), show_col_types = FALSE)
) |>
  mutate(
    sig = factor(
      case_when(
        p.value < 0.001 ~ "p < 0.001",
        p.value < 0.01  ~ "p < 0.01",
        p.value < 0.05  ~ "p < 0.05",
        TRUE            ~ "ns"
      ), levels = SIG_LEVELS
    )
  ) |>
  group_by(gene) |>
  arrange(estimate) |>
  mutate(rank = row_number()) |>
  ungroup()

sig_pct <- jk_all |>
  group_by(gene) |>
  summarise(pct = round(100 * mean(p.value < 0.05), 1), .groups = "drop")

make_jk_panel <- function(dat, gene_label, tag_letter,
                          show_x = TRUE, show_legend = FALSE, show_y = TRUE) {
  pct <- sig_pct$pct[sig_pct$gene == unique(dat$gene)]
  lbl <- paste0(gene_label, " — ", pct, "% p < 0.05")
  p <- ggplot(dat, aes(x = rank, y = estimate, colour = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(ymin = estimate - std_error,
                      ymax = estimate + std_error),
                  width = 0, alpha = 0.35, linewidth = 0.2) +
    geom_point(size = 0.6, alpha = 0.75) +
    annotate("text", x = -Inf, y = Inf, label = lbl,
             hjust = -0.05, vjust = 1.4,
             size = (BS - 1) / .pt, colour = "grey20", fontface = "bold") +
    scale_colour_manual(values = SIG_COLOURS, name = NULL,
                        guide = guide_legend(
                          override.aes = list(size = 2.5, alpha = 1),
                          nrow = 1)) +
    labs(x = if (show_x) "Sample rank (ordered by estimate)" else NULL,
         y = if (show_y) "Temperature coefficient" else NULL,
         tag = tag_letter) +
    common_theme +
    theme(
      legend.position      = if (show_legend) c(0.98, 0.05) else "none",
      legend.justification = if (show_legend) c(1, 0) else c(0, 0),
      legend.direction     = "horizontal",
      legend.background    = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x   = element_blank(),
      plot.margin          = if (show_x)
        margin(0, 5.5, 5.5, 5.5)
      else
        margin(5.5, 5.5, 0, 5.5)
    )
  if (!show_x)
    p <- p + theme(axis.text.x  = element_blank(),
                   axis.ticks.x = element_blank())
  p
}

p_jk_h1 <- make_jk_panel(filter(jk_all, gene == "H1"), "H1N1 (H1)", "E",
                          show_x = FALSE, show_legend = FALSE, show_y = FALSE)
p_jk_h3 <- make_jk_panel(filter(jk_all, gene == "H3"), "H3N2 (H3)", "",
                          show_x = TRUE,  show_legend = TRUE,  show_y = TRUE)

p_jk <- p_jk_h1 / p_jk_h3 + plot_layout(heights = c(1, 1))

# ══════════════════════════════════════════════════════════════════════════════
# Assemble
# ══════════════════════════════════════════════════════════════════════════════
row1 <- (p_map | p_scatter_h1 | p_scatter_h3) +
  plot_layout(widths = c(1.6, 1, 1))
row2 <- (p_forest | p_multilevel) + plot_layout(widths = c(1, 1))

fig_main <- row1 / row2 / p_jk +
  plot_layout(heights = c(0.6, 1.4, 1.2))

ggsave(file.path(out_dir, "fig_main.png"),
       fig_main, width = 8.27, height = 10, dpi = 200)
ggsave(file.path(out_dir, "fig_main.pdf"),
       fig_main, width = 8.27, height = 10)
message("Saved: paper/figures/fig_main.png / .pdf")
