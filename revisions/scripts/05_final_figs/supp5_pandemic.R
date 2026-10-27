# scripts/05_final_figs/supp5_pandemic.R
#
# Supplementary Figure 1: Pandemic vs post-pandemic phylolm for H1 (log_rate).
# Shows that the climate signal in the main analysis (fig1 panel C) holds when
# the pandemic years (2009–2010) and post-pandemic years (2011+) are modelled
# separately, demonstrating the result is not driven by the pandemic transition.
#
# INPUT:  data/rates/phylolm/H1_phylolm_pandemic_results.csv
# OUTPUT: paper/figures/supp5_pandemic.png / .pdf

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

BS  <- 10
BS2 <- 11

KEEP_TERMS <- c("scaled_temp", "scaled_humidity",
                "scaled_sampling_density", "scaled_delta_temp",
                "scaled_delta_humidity", "scaled_pop_density")
TERM_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative\nhumidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative\nhumidity",
  scaled_sampling_density = "Sampling\ndensity",
  scaled_pop_density      = "Population\ndensity"
)
PRED_COLOURS <- setNames(viridis(length(TERM_LABELS), end = 0.92), TERM_LABELS)

common_theme <- theme_bw(base_size = BS) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = BS + 1),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title         = element_text(size = BS2),
    axis.text          = element_text(size = BS - 1),
    legend.text        = element_text(size = BS - 2),
    legend.key.size    = unit(0.6, "lines"),
    plot.tag           = element_text(face = "bold.italic", size = BS + 2),
    plot.tag.position  = "topleft"
  )

res <- read_csv(
  file.path(PATHS$rates_phylolm, "H1_phylolm_pandemic_results.csv"),
  show_col_types = FALSE
) |>
  filter(term %in% KEEP_TERMS) |>
  mutate(
    ci_lo     = estimate - 1.96 * std_error,
    ci_hi     = estimate + 1.96 * std_error,
    sig_fdr   = !is.na(p.fdr) & p.fdr < 0.05,
    predictor = factor(TERM_LABELS[term], levels = TERM_LABELS),
    pt_fill   = ifelse(sig_fdr, PRED_COLOURS[as.character(predictor)], "white"),
    subset_label = factor(
      ifelse(subset == "pandemic", "Pandemic (2009–2010)", "Post-pandemic (2011+)"),
      levels = c("Pandemic (2009–2010)", "Post-pandemic (2011+)")
    )
  )

daic <- res |>
  distinct(subset, subset_label, AIC, AIC_null) |>
  mutate(daic_label = sprintf("ΔAIC = %.1f", AIC_null - AIC))

dodge <- position_dodge(width = 0.5)

fig <- ggplot(res, aes(x = estimate, y = subset_label,
                       colour = predictor, group = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = pt_fill),
             shape = 21, size = 2.2, stroke = 0.7, position = dodge) +
  geom_text(data = daic,
            aes(x = Inf, y = subset_label, label = daic_label),
            hjust = 1.05, vjust = -1.8, inherit.aes = FALSE,
            size = (BS - 3) / .pt, colour = "grey40") +
  scale_colour_manual(values = PRED_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev) +
  labs(x = "Regression coefficient", y = NULL) +
  common_theme +
  theme(
    legend.position      = c(0.01, 0.01),
    legend.justification = c(0, 0),
    legend.direction     = "vertical",
    legend.background    = element_rect(fill = "transparent", colour = NA),
    legend.key.spacing.y = unit(0.2, "lines"),
    legend.text          = element_text(size = BS - 3),
    plot.margin          = margin(2, 12, 2, 2, "mm")
  )

ggsave(file.path(out_dir, "supp5_pandemic.png"), fig,
       width = 130, height = 70, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp5_pandemic.pdf"), fig,
       width = 130, height = 70, units = "mm")

message("Written: supp5_pandemic.png / .pdf")