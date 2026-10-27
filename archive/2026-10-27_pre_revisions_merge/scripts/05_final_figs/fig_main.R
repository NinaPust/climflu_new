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
  library(ape)
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
                "scaled_delta_humidity", "scaled_pop_density")
TERM_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative\nhumidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative\nhumidity",
  scaled_sampling_density = "Sampling\ndensity",
  scaled_pop_density      = "Population\ndensity"
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

r2_forest <- res_main |>
  distinct(gene, subtype, protein, AIC, AIC_null) |>
  mutate(r2_label = sprintf("ΔAIC = %.1f", AIC_null - AIC))

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
  geom_text(data = r2_forest,
            aes(x = 0.05, y = protein, label = r2_label),
            hjust = 0, inherit.aes = FALSE,
            size = (BS - 3) / .pt, colour = "grey40") +
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
        legend.text          = element_text(size = BS - 3),
        legend.key.spacing.y = unit(0.1, "lines"),
        panel.grid.major.x   = element_blank(),
        aspect.ratio         = 0.5,
        strip.text           = element_blank(),
        strip.background     = element_blank()) +
  guides(colour = guide_legend(nrow = 6))

# ══════════════════════════════════════════════════════════════════════════════
# ROW 1-B: Root-to-tip regression — H1N1 and H3N2
# ══════════════════════════════════════════════════════════════════════════════
parse_clock <- function(gene) {
  path  <- file.path(PATHS$trees_timed, paste0(gene, "_treetime_out"),
                     "molecular_clock.txt")
  lines <- readLines(path)
  rate  <- as.numeric(sub(".*--rate:\\s*", "", lines[grep("rate",  lines)[1]]))
  r2    <- as.numeric(sub(".*--r\\^2:\\s*", "", lines[grep("r\\^2", lines)[1]]))
  list(rate = rate, r2 = r2)
}

make_rtt_panel <- function(gene, tag_letter, show_y = TRUE, title = NULL) {
  dates <- read_tsv(
    file.path(PATHS$trees_timed, paste0(gene, "_treetime_out"), "dates.tsv"),
    comment = "#", col_names = c("node", "date", "num_date"),
    show_col_types = FALSE
  ) |> filter(!grepl("^NODE_", node))

  clock <- parse_clock(gene)

  div_tree <- ape::read.nexus(
    file.path(PATHS$trees_timed, paste0(gene, "_treetime_out"),
              "divergence_tree.nexus")
  )
  rtt <- data.frame(
    node     = div_tree$tip.label,
    distance = ape::node.depth.edgelength(div_tree)[seq_len(ape::Ntip(div_tree))]
  )

  plot_df <- dates |>
    inner_join(rtt, by = "node") |>
    filter(is.finite(distance), is.finite(num_date))

  label <- sprintf("rate = %.2e subs/site/yr\nR² = %.2f", clock$rate, clock$r2)

  p <- ggplot(plot_df, aes(x = num_date, y = distance)) +
    geom_point(size = 0.8, alpha = 0.4, colour = viridis(1, begin = 0.5)) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.5,
                colour = "black", linetype = "dashed") +
    annotate("text", x = Inf, y = -Inf, label = label,
             hjust = 1.05, vjust = -0.3,
             size = (BS - 2) / .pt, colour = "grey30") +
    labs(
      x     = "Collection year",
      y     = if (show_y) "Root-to-tip distance (subs/site)" else NULL,
      title = title,
      tag   = tag_letter
    ) +
    common_theme +
    theme(
      aspect.ratio        = 0.8,
      plot.title          = element_text(size = BS, face = "bold",
                                         hjust = 0, colour = "grey20"),
      plot.title.position = "plot"
    )
  if (!show_y) p <- p + theme(axis.text.y  = element_blank(),
                               axis.ticks.y = element_blank())
  p
}

p_rtt_h1 <- make_rtt_panel("H1", "B", show_y = TRUE,  title = "H1N1")
p_rtt_h3 <- make_rtt_panel("H3", "",  show_y = FALSE, title = "H3N2")

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

r2_ml <- ml_combined |>
  distinct(subtype, region, AIC, AIC_null) |>
  mutate(r2_label = sprintf("ΔAIC = %.1f", AIC_null - AIC))

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
  geom_text(data = r2_ml,
            aes(x = 0.05, y = region, label = r2_label),
            hjust = 0, inherit.aes = FALSE,
            size = (BS - 3) / .pt, colour = "grey40") +
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
        aspect.ratio       = 0.65,
        strip.text         = element_blank(),
        strip.background   = element_blank())

# ══════════════════════════════════════════════════════════════════════════════
# ROWS 3 & 4: Jackknife 2×2 — whole gene (left) | epitope (right)
#             Row 1 = H1N1, Row 2 = H3N2; y-axis range matched within rows
# ══════════════════════════════════════════════════════════════════════════════
SIG_LEVELS  <- c("ns", "p < 0.05", "p < 0.01", "p < 0.001")
SIG_COLOURS <- setNames(viridis(4, end = 0.92), SIG_LEVELS)

prep_jk <- function(df) {
  df |>
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
}

jk_whole <- prep_jk(bind_rows(
  read_csv(file.path(PATHS$rates, "jackknife", "H1", "pct60", "jackknife_results.csv"),
           show_col_types = FALSE),
  read_csv(file.path(PATHS$rates, "jackknife", "H3", "pct60", "jackknife_results.csv"),
           show_col_types = FALSE)
))

jk_epi <- prep_jk(bind_rows(
  read_csv(file.path(PATHS$rates, "jackknife", "epitope",
                     "H1_epitope_jackknife_results.csv"),
           show_col_types = FALSE) |>
    rename(any_of(c(std_error = "StdErr"))) |> mutate(gene = "H1"),
  read_csv(file.path(PATHS$rates, "jackknife", "epitope",
                     "H3_epitope_jackknife_results.csv"),
           show_col_types = FALSE) |>
    rename(any_of(c(std_error = "StdErr"))) |> mutate(gene = "H3")
))

# Shared y limits per gene (whole + epitope combined)
# Both rows use the same breaks; ylims padded to ensure all breaks are visible
JK_BREAKS <- seq(-0.20, 0.05, by = 0.05)
ylims <- bind_rows(jk_whole, jk_epi) |>
  mutate(lo = estimate - std_error, hi = estimate + std_error) |>
  group_by(gene) |>
  summarise(ymin = min(lo, na.rm = TRUE), ymax = max(hi, na.rm = TRUE),
            .groups = "drop") |>
  mutate(
    ymin = pmin(ymin, min(JK_BREAKS) - 0.01),
    ymax = pmax(ymax, max(JK_BREAKS) + 0.01)
  )

make_jk_panel <- function(dat, ylim = NULL, gene_label, tag_letter,
                          show_x = TRUE, show_legend = FALSE,
                          show_y = TRUE) {
  n_draw <- unique(dat$n)
  lbl    <- paste0(gene_label, " — n = ", n_draw, " per draw")

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
    scale_y_continuous(
      limits = c(ylims$ymin[ylims$gene == unique(dat$gene)],
                 ylims$ymax[ylims$gene == unique(dat$gene)]),
      breaks = JK_BREAKS,
      labels = scales::label_number(accuracy = 0.01)
    ) +
    labs(x   = if (show_x) "Sample rank (ordered by estimate)" else NULL,
         y   = if (show_y) "Temperature coefficient" else NULL,
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
  if (!show_x) p <- p + theme(axis.text.x  = element_blank(),
                               axis.ticks.x = element_blank())
  if (!show_y) p <- p + theme(axis.text.y  = element_blank(),
                               axis.ticks.y = element_blank())
  p
}

p_jk_h1_whole <- make_jk_panel(filter(jk_whole, gene == "H1"),
                                gene_label = "H1N1", tag_letter = "E",
                                show_x = FALSE, show_y = TRUE)
p_jk_h1_epi   <- make_jk_panel(filter(jk_epi,   gene == "H1"),
                                gene_label = "H1N1 epitope", tag_letter = "",
                                show_x = FALSE, show_y = FALSE)
p_jk_h3_whole <- make_jk_panel(filter(jk_whole, gene == "H3"),
                                gene_label = "H3N2", tag_letter = "",
                                show_x = TRUE, show_y = TRUE,
                                show_legend = TRUE)
p_jk_h3_epi   <- make_jk_panel(filter(jk_epi,   gene == "H3"),
                                gene_label = "H3N2 epitope", tag_letter = "",
                                show_x = TRUE, show_y = FALSE)

p_jk <- (p_jk_h1_whole | p_jk_h1_epi) / (p_jk_h3_whole | p_jk_h3_epi) +
  plot_layout(heights = c(1, 1))

# ══════════════════════════════════════════════════════════════════════════════
# Assemble
# ══════════════════════════════════════════════════════════════════════════════
row1 <- (p_map | p_rtt_h1 | p_rtt_h3) +
  plot_layout(widths = c(1.6, 1, 1))
row2 <- (p_forest | p_multilevel) + plot_layout(widths = c(1, 1))

fig_main <- row1 / row2 / p_jk +
  plot_layout(heights = c(0.6, 1.4, 1.2))

ggsave(file.path(out_dir, "fig_main.png"),
       fig_main, width = 8.27, height = 10, dpi = 200)
ggsave(file.path(out_dir, "fig_main.pdf"),
       fig_main, width = 8.27, height = 10)

# Panel B: RTT plots occupy 2/3.6 of row-1 width, row-1 is 0.6/3.2 of total height
panel_B <- p_rtt_h1 | p_rtt_h3
ggsave(file.path(out_dir, "fig_main_B.pdf"), panel_B,
       width = 8.27 * 2 / 3.6, height = 10 * 0.6 / 3.2)

# Panel E: jackknife 2x2 occupies full width, 1.2/3.2 of total height
ggsave(file.path(out_dir, "fig_main_E.pdf"), p_jk,
       width = 8.27, height = 10 * 1.2 / 3.2)

ggsave(file.path(out_dir, "fig_main_C.pdf"), p_forest,
       width = 8.27 / 2, height = 10 * 1.4 / 3.2)

message("Saved: paper/figures/fig_main.png / .pdf / fig_main_B.pdf / fig_main_C.pdf / fig_main_E.pdf")
