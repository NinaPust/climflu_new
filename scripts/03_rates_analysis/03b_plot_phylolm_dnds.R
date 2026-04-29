# scripts/03_rates_analysis/03b_plot_phylolm_dnds.R
#
# PURPOSE: Forest plots of phylolm coefficients for per-branch dN rate,
#          dS rate, and branch omega (dN/dS).
#          Faceted by response variable and subtype (H1N1 / H3N2).
#
# INPUT:  data/rates/all_phylolm_dnds_results.csv
# OUTPUT: figures/phylolm_forest_dnds.png
#         figures/phylolm_forest_dnds_fdr.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

source(here("config.R"))

out_dir <- here("figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

PREDICTOR_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative humidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative humidity",
  scaled_sampling_density = "Sampling density",
  scaled_pop_density      = "Population density"
)

RESPONSE_LABELS <- c(
  log_dN_rate      = "dN rate",
  log_dS_rate      = "dS rate",
  log_branch_omega = "Branch ω (dN/dS)"
)

PRED_COLOURS <- setNames(
  viridisLite::viridis(length(PREDICTOR_LABELS), end = 0.92),
  PREDICTOR_LABELS
)

make_plot_df <- function(res, sig_col) {
  res |>
    filter(term %in% names(PREDICTOR_LABELS),
           response %in% names(RESPONSE_LABELS)) |>
    mutate(
      ci_lo       = estimate - 1.96 * std_error,
      ci_hi       = estimate + 1.96 * std_error,
      significant = !is.na(.data[[sig_col]]) & .data[[sig_col]] < 0.05,
      predictor   = factor(PREDICTOR_LABELS[term],
                           levels = PREDICTOR_LABELS),
      protein     = factor(
        case_when(
          gene %in% c("H1", "H3") ~ "HA",
          gene %in% c("N1", "N2") ~ "NA"
        ),
        levels = c("HA", "NA")
      ),
      subtype      = ifelse(gene %in% c("H1", "N1"), "H1N1", "H3N2"),
      response_lbl = factor(RESPONSE_LABELS[response],
                            levels = RESPONSE_LABELS),
      point_fill   = ifelse(significant,
                            PRED_COLOURS[as.character(predictor)], "white")
    ) |>
    filter(!is.na(protein))
}

make_forest <- function(plot_df, title, subtitle) {
  dodge <- position_dodge(width = 0.5)
  ggplot(plot_df,
         aes(x = estimate, y = protein,
             colour = predictor, group = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   linewidth = 0.7, position = dodge) +
    geom_point(aes(fill = point_fill),
               shape = 21, size = 3, stroke = 0.8,
               position = dodge) +
    scale_colour_manual(values = PRED_COLOURS, name = NULL) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev) +
    facet_grid(response_lbl ~ subtype) +
    labs(title = title, subtitle = subtitle,
         x = "Coefficient estimate", y = NULL) +
    theme_bw(base_size = 12) +
    theme(
      strip.text         = element_text(face = "bold", size = 11),
      legend.position    = "bottom",
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank()
    )
}

# ── Nominal p < 0.05 ──────────────────────────────────────────────────────────
res_path <- file.path(PATHS$rates, "all_phylolm_dnds_results.csv")
if (!file.exists(res_path)) stop("Results not found: ", res_path)
res <- read_csv(res_path, show_col_types = FALSE)

plot_df <- make_plot_df(res, "p.value")
p <- make_forest(
  plot_df,
  title    = "Phylogenetic regression: predictors of per-branch dN/dS rates",
  subtitle = "Coefficients ± 95% CI  |  Filled = p < 0.05"
)
out_path <- file.path(out_dir, "phylolm_forest_dnds.png")
ggsave(out_path, p, width = 10, height = 8, dpi = 180)
message("Saved: ", out_path)

# ── BH-FDR corrected ─────────────────────────────────────────────────────────
if ("p.fdr" %in% names(res)) {
  plot_df_fdr <- make_plot_df(res, "p.fdr")
  p_fdr <- make_forest(
    plot_df_fdr,
    title    = "Phylogenetic regression: predictors of per-branch dN/dS rates",
    subtitle = "Coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05"
  )
  out_path_fdr <- file.path(out_dir, "phylolm_forest_dnds_fdr.png")
  ggsave(out_path_fdr, p_fdr, width = 10, height = 8, dpi = 180)
  message("Saved: ", out_path_fdr)
}
