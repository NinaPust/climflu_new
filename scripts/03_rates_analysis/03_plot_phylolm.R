# scripts/03_rates_analysis/03_plot_phylolm.R
#
# PURPOSE: Forest plot of phylolm coefficients.
#          Panels = virus subtype (H1N1 / H3N2)
#          x-axis = protein (HA, NA, PA)
#          Colour = predictor; dodged so lines don't overlap.
#          Filled points = p < 0.05.
#
# INPUT:  data/rates/{gene}_phylolm_results.csv
# OUTPUT: figures/phylolm_forest.png
#
# USAGE:  Rscript scripts/03_rates_analysis/03_plot_phylolm.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(here("config.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

out_dir <- here("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load results ──────────────────────────────────────────────────────────────
res <- lapply(ALL_GENES, function(g) {
  path <- file.path(PATHS$rates, paste0(g, "_phylolm_results.csv"))
  if (!file.exists(path)) { warning("Missing: ", path); return(NULL) }
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

# ── Prepare plotting data ─────────────────────────────────────────────────────
PREDICTOR_LABELS <- c(
  scaled_temp        = "Temperature",
  scaled_humidity    = "Relative humidity",
  scaled_pop_density = "Population density",
  scaled_season      = "Season"
)

PRED_COLOURS <- setNames(
  viridisLite::viridis(length(PREDICTOR_LABELS), end = 0.92),
  PREDICTOR_LABELS
)

plot_df <- res |>
  filter(term %in% names(PREDICTOR_LABELS), response == "log_rate") |>
  mutate(
    ci_lo       = estimate - 1.96 * std_error,
    ci_hi       = estimate + 1.96 * std_error,
    significant = p.value < 0.05,
    predictor   = factor(PREDICTOR_LABELS[term], levels = PREDICTOR_LABELS),
    protein     = factor(
      case_when(
        gene %in% c("H1", "H3")             ~ "HA",
        gene %in% c("N1", "N2")             ~ "NA",
        gene %in% c("PA_H1N1", "PA_H3N2")  ~ "PA"
      ),
      levels = c("HA", "NA", "PA")
    ),
    subtype = ifelse(gene %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2"),
    point_fill = ifelse(significant, PRED_COLOURS[as.character(predictor)], "white")
  )

dodge <- position_dodge(width = 0.5)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(plot_df,
            aes(x = estimate, y = protein,
                colour = predictor, group = predictor)) +

  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +

  # Separator before PA (negative control)
  geom_hline(yintercept = 1.5, linetype = "dotted",
             colour = "grey70", linewidth = 0.4) +

  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 linewidth = 0.7, position = dodge) +

  geom_point(aes(fill = point_fill),
             shape = 21, size = 3, stroke = 0.8,
             position = dodge) +

  scale_colour_manual(values = PRED_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev) +

  facet_wrap(~subtype, nrow = 1) +

  labs(
    title    = "Phylogenetic regression: predictors of log terminal branch rate",
    subtitle = "Coefficients ± 95% CI  |  PA = negative control  |  Filled = p < 0.05",
    x        = "Coefficient estimate",
    y        = NULL
  ) +

  theme_bw(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", size = 13),
    legend.position    = "bottom",
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank()
  )

out_path <- file.path(out_dir, "phylolm_forest.png")
ggsave(out_path, p, width = 9, height = 6, dpi = 180)
message("Saved: ", out_path)
