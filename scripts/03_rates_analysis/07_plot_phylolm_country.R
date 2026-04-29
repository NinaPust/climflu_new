# scripts/03_rates_analysis/07_plot_phylolm_country.R
#
# PURPOSE: Forest plots for the dense-country sensitivity analysis
#          (06_phylolm_country.R). Produces two PNGs — uncorrected and
#          BH-FDR corrected — mirroring 03_plot_phylolm.R.
#
# INPUT:  data/rates/all_phylolm_country_results.csv
# OUTPUT: figures/phylolm_forest_country.png
#         figures/phylolm_forest_country_fdr.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(here("config.R"))

out_dir <- here("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

combined_path <- file.path(PATHS$rates, "all_phylolm_country_results.csv")
if (!file.exists(combined_path)) {
  stop("Combined results not found. Run 06_phylolm_country.R first:\n  ",
       combined_path)
}

res <- read_csv(combined_path, show_col_types = FALSE)

PREDICTOR_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative humidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative humidity",
  scaled_sampling_density = "Sampling density",
  scaled_pop_density      = "Population density"
  # scaled_season         = "Season"
)

PRED_COLOURS <- setNames(
  viridisLite::viridis(length(PREDICTOR_LABELS), end = 0.92),
  PREDICTOR_LABELS
)

make_plot_df <- function(res, sig_col) {
  res |>
    filter(term %in% names(PREDICTOR_LABELS), response == "log_rate") |>
    mutate(
      ci_lo       = estimate - 1.96 * std_error,
      ci_hi       = estimate + 1.96 * std_error,
      significant = !is.na(.data[[sig_col]]) & .data[[sig_col]] < 0.05,
      predictor   = factor(PREDICTOR_LABELS[term], levels = PREDICTOR_LABELS),
      protein     = factor(
        case_when(
          gene %in% c("H1", "H3")            ~ "HA",
          gene %in% c("N1", "N2")            ~ "NA",
          gene %in% c("PA_H1N1", "PA_H3N2") ~ "PA"
        ),
        levels = c("HA", "NA", "PA")
      ),
      subtype    = ifelse(gene %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2"),
      point_fill = ifelse(significant,
                          PRED_COLOURS[as.character(predictor)], "white")
    )
}

make_forest <- function(plot_df, subtitle) {
  dodge <- position_dodge(width = 0.5)

  ggplot(plot_df,
         aes(x = estimate, y = protein,
             colour = predictor, group = predictor)) +

    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +

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
      title    = paste("Phylogenetic regression (dense-country subset):",
                       "predictors of log terminal branch rate"),
      subtitle = subtitle,
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
}

# ── Uncorrected ───────────────────────────────────────────────────────────────
p_raw <- make_forest(
  make_plot_df(res, "p.value"),
  "Coefficients ± 95% CI  |  Controls: year FE + country FE  |  PA = negative control  |  Filled = p < 0.05"
)
out_raw <- file.path(out_dir, "phylolm_forest_country.png")
ggsave(out_raw, p_raw, width = 9, height = 6, dpi = 180)
message("Saved: ", out_raw)

# ── FDR-corrected ─────────────────────────────────────────────────────────────
p_fdr <- make_forest(
  make_plot_df(res, "p.fdr"),
  "Coefficients ± 95% CI  |  Controls: year FE + country FE  |  PA = negative control  |  Filled = BH-FDR p < 0.05"
)
out_fdr <- file.path(out_dir, "phylolm_forest_country_fdr.png")
ggsave(out_fdr, p_fdr, width = 9, height = 6, dpi = 180)
message("Saved: ", out_fdr)
