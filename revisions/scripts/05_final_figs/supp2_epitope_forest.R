# scripts/05_final_figs/supp2_epitope_forest.R
#
# Supplementary Figure 5: Forest plot of climate predictors of NS and S
# substitution rates across genomic regions (full gene, HA2, epitope, ridge).
# FDR-corrected significance. Layout: 2x2 grid, columns = subtype, rows = rate
# type. X-axis range is shared within columns (per subtype) but independent
# between columns.
#
# INPUT:  data/rates/phylolm/all_phylolm_ns_results.csv
# OUTPUT: paper/figures/supp2_epitope_forest.png / .pdf

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

predictor_labels <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative humidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔHumidity",
  scaled_sampling_density = "Sampling density",
  scaled_pop_density      = "Population density"
)

# REVISED SOURCE. This figure used to read all_phylolm_ns_results.csv from the
# superseded script 12, whose responses were named log_{region}_ns_rate. It now
# reads 16_phylolm_slac.R, where the region is a separate `subset` column and
# the response is one of log_dN_rate / log_dS_rate / log_omega. The old filter
# matched nothing against the new file, which surfaced as an empty data frame
# and a cryptic "limits must be a vector of length 2" error from the x scale.
RATE_TYPES <- c(log_dN_rate = "NS rate", log_dS_rate = "S rate")

REGION_LABELS_S5 <- c(full = "Full gene", HA1 = "HA1", HA2 = "HA2",
                      epitope = "Epitope", ridge = "Ridge")
REGION_LEVELS <- unname(REGION_LABELS_S5)

pred_colours <- setNames(viridis(length(predictor_labels), end = 0.92),
                         predictor_labels)

common_theme <- theme_bw(base_size = BS) +
  theme(
    panel.grid.minor     = element_blank(),
    panel.grid.major.x   = element_blank(),
    strip.background     = element_blank(),
    strip.text           = element_blank(),
    axis.title           = element_text(size = BS2),
    axis.text            = element_text(size = BS - 1),
    legend.text          = element_text(size = BS - 1),
    legend.key.size      = unit(0.6, "lines"),
    legend.key.spacing.y = unit(0.2, "lines"),
    aspect.ratio         = 1,
    plot.title           = element_text(size = BS, face = "bold", hjust = 0.5)
  )

raw <- read_csv(
  file.path(PATHS$rates_phylolm, "all_phylolm_slac_results.csv"),
  show_col_types = FALSE
) |>
  filter(response %in% names(RATE_TYPES),
         subset %in% names(REGION_LABELS_S5),
         gene %in% c("H1", "H3")) |>
  mutate(
    region    = factor(REGION_LABELS_S5[subset], levels = REGION_LEVELS),
    rate_type = unname(RATE_TYPES[response]),
    subtype   = ifelse(gene == "H1", "H1N1", "H3N2")
  )

if (nrow(raw) == 0)
  stop("No rows matched in all_phylolm_slac_results.csv — ",
       "run 04_epitope_analysis/16_phylolm_slac.R first.")

daic_df <- raw |>
  distinct(subtype, rate_type, region, AIC, AIC_null) |>
  mutate(daic_label = sprintf("ΔAIC = %.1f", AIC_null - AIC))

plot_df <- raw |>
  filter(term %in% names(predictor_labels)) |>
  mutate(
    ci_lo      = estimate - 1.96 * std_error,
    ci_hi      = estimate + 1.96 * std_error,
    sig_fdr    = !is.na(p.fdr) & p.fdr < 0.05,
    predictor  = factor(predictor_labels[term], levels = predictor_labels),
    point_fill = ifelse(sig_fdr, pred_colours[as.character(predictor)], "white")
  )

# Shared x limits per subtype (across both rate types in that column)
x_limits <- plot_df |>
  group_by(subtype) |>
  summarise(xmin = min(ci_lo, na.rm = TRUE),
            xmax = max(ci_hi, na.rm = TRUE),
            .groups = "drop")

dodge <- position_dodge(width = 0.7)

make_panel <- function(data, daic, subtype_lim, show_y, show_title,
                       col_title = NULL, show_legend = FALSE,
                       show_x = TRUE) {
  xlim <- c(subtype_lim$xmin, subtype_lim$xmax)

  ggplot(data, aes(x = estimate, y = region,
                   colour = predictor, group = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   linewidth = 0.6, position = dodge) +
    geom_point(aes(fill = point_fill), shape = 21, size = 2, stroke = 0.7,
               position = dodge) +
    geom_text(data = daic,
              aes(x = Inf, y = region, label = daic_label),
              hjust = 1.05, inherit.aes = FALSE,
              size = (BS - 3) / .pt, colour = "grey40") +
    scale_colour_manual(values = pred_colours, name = NULL) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev) +
    scale_x_continuous(limits = xlim) +
    labs(x = "Regression coefficient", y = NULL,
         title = if (show_title) col_title else NULL) +
    common_theme +
    theme(
      legend.position      = if (show_legend) c(0.01, 0.01) else "none",
      legend.justification = c(0, 0),
      legend.direction     = "vertical",
      legend.background    = element_rect(fill = "transparent", colour = NA),
      axis.text.y          = if (show_y) element_text(size = BS - 1)
                             else element_blank(),
      axis.ticks.y         = if (show_y) element_line() else element_blank(),
      axis.title.x         = if (show_x) element_text() else element_blank(),
      axis.text.x          = if (show_x) element_text() else element_blank(),
      axis.ticks.x         = if (show_x) element_line() else element_blank()
    )
}

lim_H1N1 <- filter(x_limits, subtype == "H1N1")
lim_H3N2 <- filter(x_limits, subtype == "H3N2")

p_ns_H1 <- make_panel(filter(plot_df, subtype == "H1N1", rate_type == "NS rate"),
                       filter(daic_df, subtype == "H1N1", rate_type == "NS rate"),
                       lim_H1N1, show_y = TRUE,  show_title = TRUE,
                       col_title = "H1N1", show_legend = TRUE, show_x = FALSE) +
  theme(legend.position      = c(0.01, 0.99),
        legend.justification = c(0, 1),
        legend.text          = element_text(size = BS - 2))

p_ns_H3 <- make_panel(filter(plot_df, subtype == "H3N2", rate_type == "NS rate"),
                       filter(daic_df, subtype == "H3N2", rate_type == "NS rate"),
                       lim_H3N2, show_y = FALSE, show_title = TRUE,
                       col_title = "H3N2", show_legend = FALSE, show_x = FALSE)

p_s_H1  <- make_panel(filter(plot_df, subtype == "H1N1", rate_type == "S rate"),
                       filter(daic_df, subtype == "H1N1", rate_type == "S rate"),
                       lim_H1N1, show_y = TRUE,  show_title = FALSE,
                       show_legend = FALSE)

p_s_H3  <- make_panel(filter(plot_df, subtype == "H3N2", rate_type == "S rate"),
                       filter(daic_df, subtype == "H3N2", rate_type == "S rate"),
                       lim_H3N2, show_y = FALSE, show_title = FALSE,
                       show_legend = FALSE)

# Row labels via patchwork tags on left panels
p_ns_H1 <- p_ns_H1 + labs(tag = "NS rate") +
  theme(plot.tag          = element_text(size = BS, angle = 90, hjust = 0.5,
                                         vjust = 0),
        plot.tag.position = c(0, 0.5))

p_s_H1 <- p_s_H1 + labs(tag = "S rate") +
  theme(plot.tag          = element_text(size = BS, angle = 90, hjust = 0.5,
                                         vjust = 0),
        plot.tag.position = c(0, 0.5))

fig <- (p_ns_H1 | p_ns_H3) / (p_s_H1 | p_s_H3) +
  plot_layout(guides = "keep")

ggsave(file.path(out_dir, "supp2_epitope_forest.png"), fig,
       width = 170, height = 160, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp2_epitope_forest.pdf"), fig,
       width = 170, height = 160, units = "mm")

message("Written: supp2_epitope_forest.png / .pdf")