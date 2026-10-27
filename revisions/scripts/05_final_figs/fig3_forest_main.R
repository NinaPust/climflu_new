# scripts/05_final_figs/fig3_forest_main.R
#
# Fig 3: Phylogenetic regression — main forest plot
#   Coefficients for temperature and relative humidity on log(branch rate).
#   All six genes shown (HA, NA, PA) in each subtype panel.
#   PA rows have a shaded background to highlight the negative control.
#   Filled points = BH-FDR p < 0.05.
#
# INPUT:  data/rates/phylolm/all_phylolm_results.csv
# OUTPUT: paper/figures/fig3_forest_main.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
})

source(here("revisions", "config.R"))

# Revised figures stay under revisions/paper/ so the submitted ones in
# paper/figures/ are never overwritten.
out_dir <- here("revisions", "paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load FDR results ──────────────────────────────────────────────────────────
res <- read_csv(file.path(PATHS$rates_phylolm, "all_phylolm_results.csv"),
                show_col_types = FALSE)

# ── Prepare data ──────────────────────────────────────────────────────────────
KEEP_TERMS <- c("scaled_temp", "scaled_humidity",
                "scaled_sampling_density", "scaled_delta_temp",
                "scaled_delta_humidity")

TERM_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative humidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative humidity",
  scaled_sampling_density = "Sampling density"
)

TERM_COLOURS <- setNames(
  viridis(length(TERM_LABELS), end = 0.92),
  TERM_LABELS
)

plot_df <- res |>
  filter(response == "log_rate", term %in% KEEP_TERMS) |>
  mutate(
    ci_lo      = estimate - 1.96 * std_error,
    ci_hi      = estimate + 1.96 * std_error,
    sig_fdr    = !is.na(p.fdr) & p.fdr < 0.05,
    sig_nom    = p.value < 0.05,
    protein    = factor(
      case_when(
        gene %in% c("H1", "H3")           ~ "HA",
        gene %in% c("N1", "N2")           ~ "NA",
        gene %in% c("PA_H1N1", "PA_H3N2") ~ "PA"
      ),
      levels = c("HA", "NA", "PA")
    ),
    subtype    = ifelse(gene %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2"),
    predictor  = factor(TERM_LABELS[term], levels = TERM_LABELS),
    point_fill = ifelse(sig_fdr,
                        TERM_COLOURS[as.character(predictor)],
                        "white")
  )

dodge <- position_dodge(width = 0.55)

# ── Shading rectangle for PA rows ─────────────────────────────────────────────
# With scale_y_discrete(limits = rev), factor levels c("HA","NA","PA") are
# reversed to c("PA","NA","HA"), so PA sits at y=1, NA at y=2, HA at y=3.
shade_df <- data.frame(
  subtype = c("H1N1", "H3N2"),
  ymin    = 0.5,
  ymax    = 1.5
)

p <- ggplot(plot_df,
            aes(x = estimate, y = protein,
                colour = predictor, group = predictor)) +

  # PA background shading
  geom_rect(data = shade_df,
            aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
            inherit.aes = FALSE,
            fill = "grey94", colour = NA, alpha = 1) +

  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +

  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 linewidth = 0.7, position = dodge) +

  geom_point(aes(fill = point_fill),
             shape = 21, size = 3, stroke = 0.9,
             position = dodge) +

  scale_colour_manual(values = TERM_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev) +

  facet_wrap(~ subtype, nrow = 1) +


  labs(x = "Coefficient", y = NULL) +

  theme_bw(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", size = 14),
    strip.background   = element_blank(),
    legend.position    = "bottom",
    legend.key.size    = unit(0.9, "lines"),
    legend.text        = element_text(size = 12),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    aspect.ratio       = 1,
    axis.title.x       = element_text(size = 14),
    axis.text          = element_text(size = 12),
    plot.caption       = element_text(size = 8, colour = "grey40")
  )

ggsave(file.path(out_dir, "fig3_forest_main.png"),
       p, width = 10, height = 5.5, dpi = 200)
message("Saved: paper/figures/fig3_forest_main.png")
