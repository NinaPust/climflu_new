# scripts/05_final_figs/fig5_ns_vs_s.R
#
# Fig 5: Non-synonymous and synonymous rates are both cold-associated
#   Side-by-side temperature coefficients for NS and S rates at the whole-gene
#   level, for H1 and H3. Both being negative supports a mutational (polymerase
#   fidelity / replication rate) mechanism rather than purely immune selection.
#
#   Layout:
#     x-axis: temperature coefficient (scaled)
#     y-axis: gene (H1 / H3)
#     Panels: NS rate | S rate
#     Also shows epitope NS rate for comparison in the NS panel
#
# INPUT:  data/rates/phylolm/all_phylolm_ns_results.csv
# OUTPUT: paper/figures/fig5_ns_vs_s.png

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

# ── Load ──────────────────────────────────────────────────────────────────────
ns_res <- read_csv(file.path(PATHS$rates_phylolm, "all_phylolm_ns_results.csv"),
                   show_col_types = FALSE) |>
  filter(term == "scaled_temp", gene %in% c("H1", "H3"))

# ── Shared aesthetics ─────────────────────────────────────────────────────────
GENE_COLOURS <- setNames(viridis(2, end = 0.75), c("H1", "H3"))

prep <- function(df, responses, labels) {
  df |>
    filter(response %in% responses) |>
    mutate(
      region    = factor(labels[response], levels = labels),
      ci_lo     = estimate - 1.96 * std_error,
      ci_hi     = estimate + 1.96 * std_error,
      sig_fdr   = !is.na(p.fdr) & p.fdr < 0.05,
      sig_nom   = p.value < 0.05,
      gene_col  = GENE_COLOURS[gene],
      point_fill = ifelse(sig_fdr, gene_col, "white")
    )
}

dodge <- position_dodge(width = 0.5)

make_panel <- function(df, title, xlims = NULL) {
  p <- ggplot(df, aes(x = estimate, y = region,
                       colour = gene, group = gene)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   linewidth = 1.1, position = dodge) +
    geom_point(aes(fill = point_fill),
               shape = 21, size = 4, stroke = 1,
               position = dodge) +
    scale_colour_manual(values = GENE_COLOURS, name = "Gene",
                        labels = c(H1 = "H1 (H1N1)", H3 = "H3 (H3N2)")) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev) +
    labs(x = "Coefficient", y = NULL, title = title) +
    theme_bw(base_size = 12) +
    theme(
      legend.position    = "bottom",
      legend.text        = element_text(size = 12),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      axis.title.x       = element_text(size = 14),
      axis.text          = element_text(size = 12),
      plot.title         = element_text(face = "bold", size = 14),
      aspect.ratio       = 1
    )
  if (!is.null(xlims)) p <- p + coord_cartesian(xlim = xlims)
  p
}

# NS rates: ridge, epitope, whole gene
ns_responses <- c(
  log_ridge_ns_rate   = "Ridge",
  log_epitope_ns_rate = "Epitope",
  log_full_ns_rate    = "Whole gene"
)
df_ns <- prep(ns_res, names(ns_responses), ns_responses)

# S rates: whole gene (ridge and epitope have very small n)
s_responses <- c(
  log_epitope_s_rate = "Epitope",
  log_full_s_rate    = "Whole gene"
)
df_s <- prep(ns_res, names(s_responses), s_responses)

# Align x-axis limits across panels
xlim_shared <- range(
  c(df_ns$ci_lo, df_ns$ci_hi, df_s$ci_lo, df_s$ci_hi),
  na.rm = TRUE
) * 1.05

pns <- make_panel(df_ns, "Non-synonymous (NS) rate", xlim_shared)
ps  <- make_panel(df_s,  "Synonymous (S) rate",      xlim_shared)

fig5 <- (pns | ps) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "fig5_ns_vs_s.png"),
       fig5, width = 10, height = 5, dpi = 200)
message("Saved: paper/figures/fig5_ns_vs_s.png")
