# scripts/05_final_figs/fig2_mutations.R
#
# Fig 2: Recurrent mutations and their climate associations
#   A – H1N1 ridge:   bar (independent occurrences) | climate effect heatmap
#   B – H3N2 ridge:   bar | heatmap
#   C – H1N1 epitope: bar (top 20) | heatmap
#   D – H3N2 epitope: bar (top 20) | heatmap
#
#   Y-axis order in each heatmap matches the adjacent bar plot (sorted by
#   count descending). No results survive FDR correction, showing climate
#   does not selectively drive specific amino acid changes.
#
# REVISED SOURCES. Mutation calls now come from TreeTime maximum-likelihood
# ancestral sequences (06_ml_mutations.R) rather than phangorn parsimony, and
# the climate tests from 07_wilcoxon_ml.R. The statistical design in 19 is
# unchanged from the submitted archived script 11_wilcoxon_climate.R, so any difference in this figure is
# attributable to the ancestral reconstruction rather than to the test.
#
# INPUT:  data/epitope/{H1,H3}_{ridge,epitope}_mutations_ml.csv
#         data/rates/{ridge,epitope}_mutation_climate_results_ml.csv
# OUTPUT: paper/figures/fig2_mutations.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(forcats)
  library(viridis)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS          <- 12
BS2         <- 13
TOP_N_EPITOPE <- 20
PRED_LEVELS <- c("Temperature", "Relative\nhumidity", "ΔTemperature", "ΔRelative\nhumidity")

common_theme <- theme_bw(base_size = BS) +
  theme(
    strip.background  = element_blank(),
    strip.text        = element_text(face = "bold", size = BS + 1),
    panel.grid.minor  = element_blank(),
    axis.title        = element_text(size = BS2),
    axis.text         = element_text(size = BS - 1),
    plot.tag          = element_text(face = "bold.italic", size = BS + 2),
    plot.tag.position = "topleft"
  )

HEAT_LOW  <- viridis(3, option = "viridis")[1]
HEAT_HIGH <- viridis(3, option = "plasma")[3]

# ── Data loading helpers ──────────────────────────────────────────────────────
load_counts <- function(gene, region, top_n = NULL) {
  path <- file.path(here("data", "epitope"),
                    paste0(gene, "_", region, "_mutations_ml.csv"))
  df <- read_csv(path, show_col_types = FALSE) |>
    count(change, name = "n") |>
    filter(n > 5) |>
    arrange(desc(n))
  if (!is.null(top_n)) df <- slice_head(df, n = top_n)
  # factor ordered by count ascending so highest bar is at top
  df |> mutate(change = factor(change, levels = rev(change)))
}

load_heat <- function(region, gene, mut_levels) {
  subtype_val <- ifelse(gene == "H1", "H1N1", "H3N2")
  read_csv(file.path(PATHS$rates,
                     paste0(region, "_mutation_climate_results_ml.csv")),
           show_col_types = FALSE) |>
    filter(subtype == subtype_val,
           mutation %in% mut_levels) |>
    mutate(
      sig             = !is.na(p.fdr) & p.fdr < 0.05,
      predictor_label = recode(predictor_label,
                               "Relative humidity" = "Relative\nhumidity",
                               "ΔHumidity"         = "ΔRelative\nhumidity"),
      predictor_label = factor(predictor_label, levels = PRED_LEVELS),
      mutation        = factor(mutation, levels = mut_levels)
    )
}

# ── Panel builders ────────────────────────────────────────────────────────────
make_bar <- function(counts_df, tag_letter, subtitle_text, show_y_title = TRUE) {
  ggplot(counts_df, aes(x = n, y = change)) +
    geom_col(fill = viridis(3, option = "viridis")[1], width = 0.7) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.08))) +
    labs(x = NULL,
         y = if (show_y_title) "Amino acid substitution" else NULL,
         tag = tag_letter, subtitle = subtitle_text) +
    common_theme +
    theme(panel.grid.major.y = element_blank(),
          plot.subtitle       = element_text(size = BS, colour = "grey40",
                                             face = "italic"))
}

make_heat <- function(heat_df, mut_levels, show_x_labels = TRUE) {
  ggplot(heat_df,
         aes(x = predictor_label,
             y = factor(mutation, levels = mut_levels),
             fill = effect_r)) +
    geom_tile(colour = "white", linewidth = 0.35) +
    geom_tile(data = filter(heat_df, sig),
              colour = "black", linewidth = 0.7, fill = NA) +
    geom_text(aes(label = sprintf("%.2f", effect_r)),
              size = (BS - 3) / .pt, colour = "grey20") +
    scale_fill_gradient2(
      low      = HEAT_LOW,
      mid      = "grey92",
      high     = HEAT_HIGH,
      midpoint = 0,
      limits   = c(-1, 1),
      name     = "Rank-biserial\ncorrelation",
      guide    = "none"
    ) +
    labs(x = NULL, y = NULL) +
    common_theme +
    theme(axis.text.x  = if (show_x_labels) element_text(angle = 0, hjust = 0.5)
                         else element_blank(),
          axis.ticks.x = if (show_x_labels) element_line() else element_blank(),
          axis.text.y  = element_blank(),
          axis.ticks.y = element_blank(),
          panel.grid   = element_blank(),
          legend.position = "none")
}

# ── Build each row ────────────────────────────────────────────────────────────
make_row <- function(gene, region, tag_letter, top_n = NULL,
                     show_y_title = TRUE, show_heat_x_labels = TRUE) {
  subtype_label <- ifelse(gene == "H1", "H1N1", "H3N2")
  region_label  <- ifelse(region == "ridge", "Ridge", "Epitope")
  subtitle      <- paste0(subtype_label, " – ", region_label,
                          if (!is.null(top_n))
                            paste0(" (top ", top_n, " mutations)") else "")

  counts     <- load_counts(gene, region, top_n = top_n)
  mut_levels <- levels(counts$change)
  heat       <- load_heat(region, gene, mut_levels)

  p_bar  <- make_bar(counts, tag_letter, subtitle, show_y_title = show_y_title)
  p_heat <- make_heat(heat, mut_levels, show_x_labels = show_heat_x_labels)

  (p_bar | p_heat) + plot_layout(widths = c(1, 2))
}

row_A <- make_row("H1", "ridge",   "A", show_heat_x_labels = FALSE)
row_B <- make_row("H3", "ridge",   "B", show_y_title = FALSE, show_heat_x_labels = FALSE)
row_C <- make_row("H1", "epitope", "C", top_n = TOP_N_EPITOPE)
row_D <- make_row("H3", "epitope", "D", top_n = TOP_N_EPITOPE, show_y_title = FALSE)

# ── Standalone legend ─────────────────────────────────────────────────────────
legend_plot <- ggplot(data.frame(x = 0, y = 0, z = 0),
                      aes(x, y, fill = z)) +
  geom_tile() +
  scale_fill_gradient2(
    low = HEAT_LOW, mid = "grey92", high = HEAT_HIGH,
    midpoint = 0, limits = c(-1, 1), name = "Rank-biserial\ncorrelation",
    guide = guide_colourbar(barwidth = 7, barheight = 0.5,
                            title.position = "top", title.hjust = 0.5,
                            direction = "horizontal")
  ) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = BS - 1),
        legend.title    = element_text(size = BS, hjust = 0.5))

legend_grob <- cowplot::get_legend(legend_plot)

sig_df <- data.frame(
  x     = c(1, 1),
  y     = c(2, 1),
  sig   = c(TRUE, FALSE),
  label = c("FDR p < 0.05", "n.s.")
)

sig_legend <- ggplot(sig_df, aes(x, y)) +
  geom_tile(fill = "grey88", colour = "white", linewidth = 0.35,
            width = 0.8, height = 0.8) +
  geom_tile(data = sig_df[sig_df$sig, ],
            fill = "grey88", colour = "black", linewidth = 0.8,
            width = 0.8, height = 0.8) +
  geom_text(aes(label = label), x = 1.55, hjust = 0,
            size = (BS - 3) / .pt, colour = "grey20") +
  xlim(0.5, 4) + ylim(0.5, 2.5) +
  theme_void()

legend_col <- cowplot::plot_grid(legend_grob, sig_legend,
                                 ncol = 1, rel_heights = c(1, 0.35))

# ── Assemble ──────────────────────────────────────────────────────────────────
n_ridge <- max(nrow(load_counts("H1", "ridge")),
               nrow(load_counts("H3", "ridge")))

main_grid <- (wrap_elements(row_A) | wrap_elements(row_B)) /
             (wrap_elements(row_C) | wrap_elements(row_D)) +
  plot_layout(heights = c(n_ridge, TOP_N_EPITOPE))

fig2 <- cowplot::plot_grid(main_grid, legend_col,
                           ncol = 2, rel_widths = c(1, 0.1))

ggsave(file.path(out_dir, "fig2_mutations.png"),
       fig2, width = 17, height = 11, dpi = 200, bg = "white")
ggsave(file.path(out_dir, "fig2_mutations.pdf"),
       fig2, width = 17, height = 11, bg = "white")
message("Saved: paper/figures/fig2_mutations.png / .pdf")
