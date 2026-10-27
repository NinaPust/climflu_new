# scripts/05_final_figs/supp7_epitope_jackknife.R
#
# Supplementary Figure 6: Jackknife robustness of temperature effect on NS
# substitution rates across genomic regions. Layout mirrors supp5: columns =
# gene (H1N1, H3N2), rows = region (HA2, Epitope, Ridge). Each panel shows
# jackknife samples ordered by estimate, coloured by significance. Full gene
# jackknife is shown separately in supp3.
#
# REVISED SOURCE: 04_epitope_analysis/21_jackknife_slac.R, which fits the full
# model (six climate predictors + year) on SLAC-derived dN and stores every
# coefficient. The submitted version read script 13's output, which jackknifed a
# temperature-only formula on parsimony NS rates.
#
# INPUT:
#   revisions/data/rates/jackknife/slac/{draw}/{gene}_{region}_{response}_results.csv
#
# OUTPUT: paper/figures/supp7_epitope_jackknife.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(viridis)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "jackknife_utils.R"))

# Revised figures stay under revisions/paper/ so the submitted ones in
# paper/figures/ are never overwritten.
out_dir <- here("revisions", "paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS  <- 9
BS2 <- 10

REGIONS <- c("HA2", "epitope", "ridge")
REGION_LABELS <- c(HA2 = "HA2", epitope = "Epitope", ridge = "Ridge")
GENES   <- c("H1", "H3")

SIG_LEVELS  <- c("ns", "*", "**", "***")
SIG_LABELS  <- c("ns" = "p > 0.05", "*" = "p < 0.05",
                 "**" = "p < 0.01", "***" = "p < 0.001")
SIG_COLOURS <- setNames(viridis(4, end = 0.92), SIG_LEVELS)
JK_COEF <- "scaled_temp"

classify_sig <- function(p) {
  factor(
    case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    levels = SIG_LEVELS
  )
}

common_theme <- theme_bw(base_size = BS) +
  theme(
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = BS2),
    axis.text          = element_text(size = BS - 1),
    legend.text        = element_text(size = BS - 2),
    legend.key.size    = unit(0.5, "lines"),
    legend.key.spacing.y = unit(0.1, "lines"),
    plot.title         = element_text(size = BS, face = "bold", hjust = 0.5),
    aspect.ratio       = 0.35
  )

# Draw size and response to plot; must match a run of 21_jackknife_slac.R.
JK_DRAW     <- "pct60"
JK_RESPONSE <- "log_dN_rate"
JK_COEF     <- "scaled_temp"

jk_dir <- file.path(PATHS$rates, "jackknife", "slac", JK_DRAW)

load_results <- function(gene, region) {
  path <- file.path(jk_dir,
                    paste0(gene, "_", region, "_", JK_RESPONSE, "_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_jackknife(path, JK_COEF, quiet = TRUE) |>
    mutate(sig = classify_sig(p.value)) |>
    arrange(estimate) |>
    mutate(rank = row_number(), gene = gene, region = region)
}

all_data <- lapply(REGIONS, function(reg) {
  lapply(GENES, function(g) load_results(g, reg)) |> bind_rows()
}) |> bind_rows()

# Only keep regions that script 21 has actually been run for; a missing region
# would otherwise reach min() on an empty vector and fail inside the x scale.
have <- all_data |> count(region) |> pull(region)
missing <- setdiff(REGIONS, have)
if (length(missing))
  message("Not yet run in 21_jackknife_slac.R (", JK_DRAW, "/", JK_RESPONSE,
          "), skipping: ", paste(missing, collapse = ", "),
          "\n  -> Rscript revisions/scripts/04_epitope_analysis/21_jackknife_slac.R",
          " --region ", paste(missing, collapse = " "))
REGIONS <- intersect(REGIONS, have)
if (length(REGIONS) == 0)
  stop("No jackknife results found under ", jk_dir,
       " — run 04_epitope_analysis/21_jackknife_slac.R first.")

x_limits <- all_data |>
  mutate(ci_lo = estimate - std_error, ci_hi = estimate + std_error) |>
  group_by(region) |>
  summarise(xmin = min(ci_lo, na.rm = TRUE),
            xmax = max(ci_hi, na.rm = TRUE),
            .groups = "drop")

make_panel <- function(gene, region, show_x = FALSE, show_y = FALSE,
                       show_legend = FALSE, col_title = NULL,
                       row_label = NULL) {
  dat   <- filter(all_data, gene == !!gene, region == !!region) |>
    mutate(ci_lo = estimate - std_error, ci_hi = estimate + std_error)
  xlim  <- x_limits |> filter(region == !!region)
  n_draw <- unique(dat$n)

  p <- ggplot(dat, aes(x = rank, y = estimate, colour = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                  width = 0.4, alpha = 0.5, linewidth = 0.3) +
    geom_point(size = 0.8, alpha = 0.8) +
    annotate("text", x = -Inf, y = Inf,
             label = paste0("n = ", n_draw, " per draw"),
             hjust = -0.05, vjust = 1.4,
             size = (BS - 1) / .pt, colour = "grey40") +
    scale_colour_manual(name = NULL, values = SIG_COLOURS,
                        labels = SIG_LABELS) +
    scale_y_continuous(limits = c(xlim$xmin, xlim$xmax)) +
    labs(
      x     = if (show_x) "Sample (ordered by estimate)" else NULL,
      y     = NULL,
      title = col_title
    ) +
    common_theme +
    theme(
      axis.text.x  = if (show_x) element_text() else element_blank(),
      axis.ticks.x = if (show_x) element_line() else element_blank(),
      axis.text.y  = if (show_y) element_text() else element_blank(),
      axis.ticks.y = if (show_y) element_line() else element_blank(),
      legend.position      = if (show_legend) c(0.99, 0.01) else "none",
      legend.justification = c(1, 0),
      legend.direction     = "horizontal",
      legend.background    = element_rect(fill = "transparent", colour = NA)
    )
  p
}

# Build 3 × 2 grid: rows = HA2 / Epitope / Ridge, columns = H1N1 / H3N2
panels <- list()
for (i in seq_along(REGIONS)) {
  reg      <- REGIONS[i]
  is_top   <- i == 1
  is_bot   <- i == length(REGIONS)

  p_H1 <- make_panel("H1", reg,
                      show_x      = is_bot,
                      show_y      = TRUE,
                      show_legend = is_top,
                      col_title   = if (is_top) "H1N1" else NULL)

  p_H3 <- make_panel("H3", reg,
                      show_x    = is_bot,
                      show_y    = FALSE,
                      col_title = if (is_top) "H3N2" else NULL)

  # Row label as rotated tag on the right side of the H3 panel
  p_H3 <- p_H3 + labs(tag = REGION_LABELS[reg]) +
    theme(plot.tag          = element_text(size = BS, angle = 270,
                                           hjust = 0.5, vjust = 0),
          plot.tag.position = c(1, 0.5))

  panels[[paste(reg, "H1", sep = "_")]] <- p_H1
  panels[[paste(reg, "H3", sep = "_")]] <- p_H3
}

# Rows are assembled from whatever regions are present rather than hardcoded,
# so the figure still builds when 21_jackknife_slac.R has only been run for a
# subset of regions.
rows <- lapply(REGIONS, function(reg)
  panels[[paste(reg, "H1", sep = "_")]] | panels[[paste(reg, "H3", sep = "_")]])
fig <- Reduce(`/`, rows) + plot_layout(guides = "keep")

# Shared y-axis title
fig_labelled <- wrap_elements(full = fig) +
  labs(tag = "Regression coefficient of temperature") +
  theme(
    plot.tag          = element_text(size = BS, angle = 90,
                                     hjust = 0.5, vjust = 1),
    plot.tag.position = c(0.03, 0.5)
  )

fig_h <- 55 * length(REGIONS)   # mm per region row
ggsave(file.path(out_dir, "supp7_epitope_jackknife.png"), fig_labelled,
       width = 290, height = fig_h, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp7_epitope_jackknife.pdf"), fig_labelled,
       width = 290, height = fig_h, units = "mm")

message("Written: supp7_epitope_jackknife.png / .pdf")