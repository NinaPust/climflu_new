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

source(here("revisions", "config.R"))
source(here("revisions", "R", "jackknife_utils.R"))

# Revised figures go under revisions/paper/ so the submitted figures in
# paper/figures/ are left untouched and the two can be compared side by side.
out_dir <- here("revisions", "paper", "figures")
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
# Sampling locations are unchanged by the revision, so these come from the
# ORIGINAL data tree. revisions/data/rates/ holds only the re-derived tables
# ({gene}_rates_slac.csv etc.), not the base {gene}_rates.csv.
rates_all <- lapply(ALL_GENES, function(g) {
  p <- here("data", "rates", paste0(g, "_rates.csv"))
  if (!file.exists(p)) return(NULL)
  read_csv(p, show_col_types = FALSE) |> mutate(gene = g)
}) |> bind_rows()
if (nrow(rates_all) == 0)
  stop("No rates tables found under data/rates/ — panel A needs them for the map")

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

# Panel C is the whole-gene phylolm from 03_rates_analysis/02_phylolm.R, which
# the revision did not touch — read from the ORIGINAL data tree.
res_main <- read_csv(
  here("data", "rates", "phylolm", "all_phylolm_results.csv"),
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
# PANEL D: Multi-level forest — temperature coefficient
#
# Layout is unchanged from the submitted figure: two facets (H1N1, H3N2), six
# rows each, narrowing from the whole gene down to the antigenic ridge:
#
#   Terminal branch rate   whole gene, overall substitutions  <- same source as C
#   Whole-gene dN          whole gene, non-synonymous
#   Whole-gene dS          whole gene, synonymous
#   HA2 dN                 conserved subunit
#   Epitope dN             antigenic sites
#   Ridge dN               antigenic ridge
#
# Only the SOURCES change. Row 1 still comes from 03_rates_analysis/02_phylolm.R
# (unrevised, read from the original data tree). Rows 2-6 now come from
# 04_epitope_analysis/16_phylolm_slac.R, i.e. HyPhy SLAC maximum-likelihood
# counts normalised by per-branch expected site counts, in place of the
# parsimony NS/S rates the old script 12 produced.
# ══════════════════════════════════════════════════════════════════════════════
ml_main <- read_csv(
  here("data", "rates", "phylolm", "all_phylolm_results.csv"),
  show_col_types = FALSE
) |>
  filter(response == "log_rate", term == "scaled_temp", gene %in% c("H1", "H3")) |>
  mutate(region = "Terminal\nbranch rate")

# Rows 2-6: one SLAC model each, selected by (region, response) pair.
SLAC_ROWS <- tibble::tribble(
  ~subset,   ~response,      ~region,
  "full",    "log_dN_rate",  "Whole-gene\ndN",
  "full",    "log_dS_rate",  "Whole-gene\ndS",
  "HA2",     "log_dN_rate",  "HA2\ndN",
  "epitope", "log_dN_rate",  "Epitope\ndN",
  "ridge",   "log_dN_rate",  "Ridge\ndN"
)

ml_slac <- read_csv(
  file.path(PATHS$rates_phylolm, "all_phylolm_slac_results.csv"),
  show_col_types = FALSE
) |>
  filter(term == "scaled_temp", gene %in% c("H1", "H3")) |>
  inner_join(SLAC_ROWS, by = c("subset", "response"))

ROW_LEVELS <- c("Terminal\nbranch rate", SLAC_ROWS$region)

ml_combined <- bind_rows(ml_main, ml_slac) |>
  mutate(
    subtype = ifelse(gene == "H1", "H1N1", "H3N2"),
    region  = factor(region, levels = ROW_LEVELS),
    ci_lo   = estimate - 1.96 * std_error,
    ci_hi   = estimate + 1.96 * std_error,
    sig_fdr = !is.na(p.fdr) & p.fdr < 0.05,
    rate_type = factor(case_when(
      region == "Terminal\nbranch rate" ~ "Substitution\nrate",
      grepl("dS$", as.character(region)) ~ "Synonymous\nrate",
      TRUE                               ~ "Non-synonymous\nrate"
    ), levels = c("Substitution\nrate", "Non-synonymous\nrate",
                  "Synonymous\nrate"))
  ) |>
  filter(!is.na(region))

TYPE_COLOURS <- setNames(viridis(3, end = 0.85),
                         c("Substitution\nrate", "Non-synonymous\nrate",
                           "Synonymous\nrate"))

# One model per row, so delta AIC is well defined per row.
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
  # Anchored to the panel edge rather than a fixed x: the SLAC coefficients
  # span a wider range than the submitted ones, so a hardcoded x lands on the
  # data for some rows.
  geom_text(data = r2_ml,
            aes(x = Inf, y = region, label = r2_label),
            hjust = 1.05, inherit.aes = FALSE,
            size = (BS - 3) / .pt, colour = "grey40") +
  scale_colour_manual(values = TYPE_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev, expand = expansion(add = 0.4)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.32))) +
  facet_wrap(~ subtype, ncol = 1) +
  geom_text(data = data.frame(subtype = "H3N2"),
            aes(x = -Inf, y = -Inf,
                label = "● FDR < 0.05   ○ n.s."),
            hjust = -0.05, vjust = -0.4, inherit.aes = FALSE,
            size = (BS - 3) / .pt, colour = "grey30") +
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

# Both jackknives are now REVISED and fit the full model (six climate
# predictors + year factor), so both store every coefficient, one row per
# replicate per term. read_jackknife() reduces each to the coefficient being
# plotted — without it a 1000-replicate run reads as 6000 samples.
JK_COEF <- "scaled_temp"

# Whole gene: 03_rates_analysis/07_phylolm_jackknife.R, 60% draw size.
jk_whole <- prep_jk(bind_rows(lapply(c("H1", "H3"), function(g) {
  p <- file.path(PATHS$rates, "jackknife", g, "pct60",
                 "jackknife_results.csv")
  if (!file.exists(p))
    stop("Whole-gene jackknife not found: ", p,
         "\nRun 03_rates_analysis/07_phylolm_jackknife.R first.")
  read_jackknife(p, JK_COEF) |> mutate(gene = g)
})))

# Epitope: 04_epitope_analysis/21_jackknife_slac.R, dN response, 60% draw —
# the same draw size as the whole-gene panel above, so the two rows are
# comparable. The submitted version used a 60% draw capped at 1000, which left
# H3 subsampling at 46% while H1 was at 60%.
JK_RESPONSE <- "log_dN_rate"
JK_DRAW     <- "pct60"
jk_epi <- prep_jk(bind_rows(lapply(c("H1", "H3"), function(g) {
  p <- file.path(PATHS$rates, "jackknife", "slac", JK_DRAW,
                 paste0(g, "_epitope_", JK_RESPONSE, "_results.csv"))
  if (!file.exists(p))
    stop("Epitope jackknife not found: ", p,
         "\nRun 04_epitope_analysis/21_jackknife_slac.R first.")
  read_jackknife(p, JK_COEF) |> mutate(gene = g)
})))

# Y range is shared WITHIN each row, not across the whole panel. The whole-gene
# and epitope coefficients differ by roughly an order of magnitude (whole gene
# spans about -0.08 to 0, epitope about -0.30 to 0), so forcing one scale on
# both flattens the whole-gene series into a line. Each row therefore gets its
# own limits, and the two genes within a row share them so they stay comparable.
#' @param top  Fixed upper limit. Both rows top out at the same value so the
#'   zero line sits at a comparable height and the four facets read as one
#'   panel, even though their lower limits differ.
#' @param grid Break spacing; NULL picks a value giving ~6 gridlines.
#' @param bottom Fixed lower limit. Given with `top`, breaks run exactly from
#'   bottom to top, so the label count is deterministic rather than depending on
#'   whether a computed break falls inside the padded range.
y_scale_for <- function(df, pad_frac = 0.03, grid = NULL,
                        top = NULL, bottom = NULL) {
  v   <- c(df$estimate - df$std_error, df$estimate + df$std_error)
  raw <- range(v, na.rm = TRUE)
  pad <- pad_frac * diff(raw)
  if (is.null(grid)) {
    cand <- c(0.005, 0.01, 0.02, 0.025, 0.05, 0.1)
    grid <- cand[which.min(abs(diff(raw) / cand - 6))]
  }
  hi <- if (is.null(top))    raw[2] + pad else top
  lo <- if (is.null(bottom)) raw[1] - pad else bottom
  if (!is.null(bottom) && raw[1] < bottom)
    warning("y lower limit ", bottom, " clips data reaching ", round(raw[1], 4))

  brk <- if (!is.null(top) && !is.null(bottom)) {
    seq(bottom, top, by = grid)
  } else {
    seq(floor(raw[1] / grid) * grid, ceiling(hi / grid) * grid, by = grid)
  }
  # Label precision from the grid's magnitude: a 0.03 grid still wants two
  # decimals, so accuracy must not simply be the grid value.
  list(lim = c(lo, hi), breaks = brk,
       acc = 10^floor(log10(grid)))
}

# Each row gets a round ceiling and a break interval matched to its own span:
# whole-gene coefficients run to about -0.10, epitope to about -0.37, so the
# epitope row needs the coarser grid to avoid a crowded axis.
# Five labels on each row. Whole gene: -0.10 to 0.02 at 0.03 gives
# -0.10 -0.07 -0.04 -0.01 0.02, all two-decimal, and -0.10 clears the lowest
# datum (-0.097). Epitope: 0.1 spacing lands on -0.3 -0.2 -0.1 0.0 0.1.
Y_WHOLE <- y_scale_for(jk_whole, grid = 0.03, top = 0.02, bottom = -0.10)
Y_EPI   <- y_scale_for(jk_epi,   grid = 0.1,  top = 0.1)

message("Panel E y ranges — whole gene: [",
        round(Y_WHOLE$lim[1], 3), ", ", round(Y_WHOLE$lim[2], 3),
        "]   epitope: [", round(Y_EPI$lim[1], 3), ", ",
        round(Y_EPI$lim[2], 3), "]")

# show_ylab is separate from show_y: both stacked rows need y tick labels, but
# only one should carry the axis title, otherwise patchwork draws it twice and
# the two overlap.
make_jk_panel <- function(dat, yscale, gene_label, tag_letter,
                          show_x = TRUE, show_legend = FALSE,
                          show_y = TRUE, show_ylab = show_y) {
  n_draw  <- unique(dat$n)
  pct_sig <- 100 * mean(dat$p.value < 0.05, na.rm = TRUE)
  lbl <- sprintf("%s — n = %d per draw\n%.1f%% of draws p < 0.05",
                 gene_label, n_draw, pct_sig)

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
      limits = yscale$lim,
      breaks = yscale$breaks,
      labels = scales::label_number(accuracy = yscale$acc)
    ) +
    labs(x   = if (show_x) "Sample rank (ordered by estimate)" else NULL,
         y   = if (show_ylab) "Temperature coefficient" else NULL,
         tag = tag_letter) +
    common_theme +
    theme(
      legend.position      = if (show_legend) c(0.98, 0.05) else "none",
      legend.justification = if (show_legend) c(1, 0) else c(0, 0),
      legend.direction     = "horizontal",
      legend.background    = element_rect(fill = "transparent", colour = NA),
      panel.grid.major.x   = element_blank(),
      # Each row of panel E is only ~1.9in tall; the default axis title size
      # overflows that and gets clipped.
      axis.title.y         = element_text(size = BS - 2),
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

# Layout is analysis-major: whole gene on the top row, epitope on the bottom,
# each row with its own y scale. Within a row the right-hand (H3N2) panel drops
# its y tick labels because it shares the left panel's scale.
p_jk_h1_whole <- make_jk_panel(filter(jk_whole, gene == "H1"), Y_WHOLE,
                               gene_label = "H1N1 whole gene", tag_letter = "E",
                               show_x = FALSE, show_y = TRUE)
p_jk_h3_whole <- make_jk_panel(filter(jk_whole, gene == "H3"), Y_WHOLE,
                               gene_label = "H3N2 whole gene", tag_letter = "",
                               show_x = FALSE, show_y = FALSE)
p_jk_h1_epi   <- make_jk_panel(filter(jk_epi,   gene == "H1"), Y_EPI,
                               gene_label = "H1N1 epitope", tag_letter = "",
                               show_x = TRUE, show_y = TRUE)
p_jk_h3_epi   <- make_jk_panel(filter(jk_epi,   gene == "H3"), Y_EPI,
                               gene_label = "H3N2 epitope", tag_letter = "",
                               show_x = TRUE, show_y = FALSE,
                               show_legend = TRUE)

p_jk <- (p_jk_h1_whole | p_jk_h3_whole) / (p_jk_h1_epi | p_jk_h3_epi) +
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

# Panel D: same geometry as C — the two share row 2 at half width each
ggsave(file.path(out_dir, "fig_main_D.pdf"), p_multilevel,
       width = 8.27 / 2, height = 10 * 1.4 / 3.2)

message("Saved to revisions/paper/figures/: fig_main.png / .pdf / fig_main_B.pdf / fig_main_C.pdf / fig_main_D.pdf / fig_main_E.pdf")
