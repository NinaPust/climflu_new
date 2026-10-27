# scripts/05_final_figs/supp4_jackknife_sequences.R
#
# Supplementary Figure 3: Per-sequence jackknife effect size distributions
# for H1 and H3. Layout: H1 (left column) | H3 (right column), one row per
# draw size (pct15, pct30, pct60). Each panel shows sequences ordered by
# median temperature coefficient, with IQR and 5-95th percentile ribbons
# coloured by sampling frequency (shared colour scale across all panels).
#
# INPUT:
#   data/rates/jackknife/{H1,H3}/{pct15,pct30,pct60}/jackknife_results.csv
#   data/rates/jackknife/{H1,H3}/{pct15,pct30,pct60}/jackknife_sequence_membership.csv
#
# OUTPUT: paper/figures/supp4_jackknife_sequences.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
})

source(here("config.R"))
source(here("R", "jackknife_utils.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS <- 9

draw_sizes <- list(
  list(label = "n1000", subdir = NULL),
  list(label = "pct15", subdir = "pct15"),
  list(label = "pct30", subdir = "pct30"),
  list(label = "pct60", subdir = "pct60")
)

genes <- c("H1", "H3")
JK_COEF <- "scaled_temp"

load_seq_summary <- function(gene, ds) {
  jk_dir <- file.path(PATHS$rates, "jackknife", gene)
  if (!is.null(ds$subdir)) jk_dir <- file.path(jk_dir, ds$subdir)
  res_path <- file.path(jk_dir, "jackknife_results.csv")
  mem_path <- file.path(jk_dir, "jackknife_sequence_membership.csv")
  if (!file.exists(res_path) || !file.exists(mem_path)) return(NULL)

  results <- read_jackknife(res_path, JK_COEF)
  if ("StdErr" %in% names(results) && !"std_error" %in% names(results))
    results <- rename(results, std_error = StdErr)
  membership <- read_csv(mem_path, show_col_types = FALSE)

  n_runs  <- nrow(results)
  n_draw  <- unique(results$n)[1]
  membership |>
    left_join(results |> select(sample_id, estimate), by = "sample_id") |>
    group_by(accession) |>
    summarise(
      median_est  = median(estimate, na.rm = TRUE),
      q25         = quantile(estimate, 0.25, na.rm = TRUE),
      q75         = quantile(estimate, 0.75, na.rm = TRUE),
      q05         = quantile(estimate, 0.05, na.rm = TRUE),
      q95         = quantile(estimate, 0.95, na.rm = TRUE),
      sample_freq = n() / n_runs * 100,
      .groups     = "drop"
    ) |>
    arrange(median_est) |>
    mutate(
      rank       = row_number(),
      gene       = gene,
      draw       = ds$label,
      n_draw     = n_draw,
      sample_freq = (sample_freq - min(sample_freq)) /
                    (max(sample_freq) - min(sample_freq))
    )
}

# Load all data to compute a global colour scale
all_summaries <- lapply(draw_sizes, function(ds) {
  lapply(genes, load_seq_summary, ds = ds) |> bind_rows()
}) |> bind_rows()

make_seq_panel <- function(dat, show_x = FALSE, show_legend = FALSE) {
  gene         <- unique(dat$gene)
  draw         <- unique(dat$draw)
  n_draw       <- unique(dat$n_draw)
  pct_str <- switch(draw,
    n1000 = "", pct15 = " 15%", pct30 = " 30%", pct60 = " 60%"
  )
  panel_label <- paste0(gene, pct_str, " (n = ", n_draw, ")")
  crosses_zero <- min(dat$q05, na.rm = TRUE) < 0 &
                  max(dat$q95, na.rm = TRUE) > 0

  p <- ggplot(dat, aes(x = rank, fill = sample_freq, colour = sample_freq))
  if (crosses_zero)
    p <- p + geom_hline(yintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4)
  p +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25, colour = NA) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.45, colour = NA) +
    geom_line(aes(y = median_est), linewidth = 0.5) +
    scale_fill_viridis_c(
      name   = "Relative sampling frequency",
      option = "D",
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("Low", "High"),
      guide  = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(0.3, "cm"))
    ) +
    scale_colour_viridis_c(
      name   = "Relative sampling frequency",
      option = "D",
      limits = c(0, 1),
      breaks = c(0, 1),
      labels = c("Low", "High"),
      guide  = guide_colourbar(barwidth = unit(3, "cm"), barheight = unit(0.3, "cm"))
    ) +
    annotate("text", x = -Inf, y = Inf,
             label = panel_label,
             hjust = -0.05, vjust = 1.4,
             size = 9.5 / .pt, colour = "black") +
    labs(x = NULL, y = NULL) +
    theme_bw(base_size = BS) +
    theme(
      panel.grid.minor      = element_blank(),
      axis.text             = element_text(size = 9),
      axis.title            = element_text(size = BS + 1),
      axis.text.x           = if (show_x) element_text(size = 9) else element_blank(),
      axis.ticks.x          = if (show_x) element_line() else element_blank(),
      plot.margin           = margin(2, 4, 1, 4, "mm"),
      legend.position       = if (show_legend) c(0.06, 0.01) else "none",
      legend.justification  = c(0, 0),
      legend.direction      = "horizontal",
      legend.background     = element_rect(fill = "transparent", colour = NA),
      legend.title          = element_text(size = 9),
      legend.title.position = "top",
      legend.text           = element_text(size = 9),
      legend.key.size       = unit(0.45, "lines")
    )
}

# Build column of panels for one gene; x-axis ticks/labels on the last row.
# show_legend_last: show the legend inside the last panel.
build_gene_column <- function(gene, show_legend_last = FALSE) {
  panels <- list()

  for (ds in draw_sizes) {
    dat <- all_summaries |>
      filter(.data$gene == .env$gene, draw == ds$label)
    if (nrow(dat) == 0) {
      message("  [", gene, " ", ds$label, "] not found, skipping")
      next
    }
    panels[[ds$label]] <- make_seq_panel(dat, show_x = FALSE)
  }

  if (length(panels) == 0) stop("No data found for ", gene)

  # Rebuild last panel with x-axis ticks and optionally the legend
  last     <- tail(names(panels), 1)
  last_dat <- all_summaries |> filter(.data$gene == .env$gene, draw == last)
  panels[[last]] <- make_seq_panel(last_dat, show_x = TRUE,
                                   show_legend = show_legend_last)

  wrap_plots(panels, ncol = 1)
}

message("Building H1 column...")
col_h1 <- build_gene_column("H1", show_legend_last = TRUE)

message("Building H3 column...")
col_h3 <- build_gene_column("H3", show_legend_last = FALSE)

fig <- (col_h1 | col_h3) +
  plot_annotation(theme = theme(
    plot.margin = margin(t = 5, r = 10, b = 15, l = 20, unit = "mm")
  ))

ggsave(file.path(out_dir, "supp4_jackknife_sequences.png"), fig,
       width = 170, height = 250, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp4_jackknife_sequences.pdf"), fig,
       width = 170, height = 250, units = "mm")

message("Written: supp4_jackknife_sequences.png / .pdf")
