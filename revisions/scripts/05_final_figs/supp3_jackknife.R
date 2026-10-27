# scripts/05_final_figs/supp3_jackknife.R
#
# Supplementary Figure 8: Jackknife subsampling sensitivity — H1 and H3
# side by side, all draw sizes stacked, shared legend.
#
# Panels in the same row (same draw size) share a common y-axis range.
# The legend appears once, in the bottom-right corner of the first panel
# (H1 n1000), styled to match panel E of the main figure.
#
# INPUT:
#   data/rates/jackknife/{gene}/jackknife_results.csv
#   data/rates/jackknife/{gene}/{pct*}/jackknife_results.csv
#
# OUTPUT:
#   paper/figures/supp3_jackknife.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "jackknife_utils.R"))

# Revised figures stay under revisions/paper/ so the submitted ones in
# paper/figures/ are never overwritten.
out_dir <- here("revisions", "paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

draw_sizes <- list(
  list(label = "n1000", subdir = NULL),
  list(label = "pct15", subdir = "pct15"),
  list(label = "pct30", subdir = "pct30"),
  list(label = "pct60", subdir = "pct60")
)

genes <- c("H1", "H3")

sig_levels  <- c("ns", "p < 0.05", "p < 0.01", "p < 0.001")
sig_colours <- setNames(viridis(4, end = 0.92), sig_levels)

classify_sig <- function(p) {
  factor(
    case_when(
      p < 0.001 ~ "p < 0.001",
      p < 0.01  ~ "p < 0.01",
      p < 0.05  ~ "p < 0.05",
      TRUE      ~ "ns"
    ),
    levels = sig_levels
  )
}

# ── Pre-load all results and compute per-draw-size y limits ──────────────────
load_results <- function(gene, ds) {
  jk_dir <- file.path(PATHS$rates, "jackknife", gene)
  if (!is.null(ds$subdir)) jk_dir <- file.path(jk_dir, ds$subdir)
  res_path <- file.path(jk_dir, "jackknife_results.csv")
  if (!file.exists(res_path)) return(NULL)
  res <- read_jackknife(res_path, JK_COEF)
  if ("StdErr" %in% names(res) && !"std_error" %in% names(res))
    res <- rename(res, std_error = StdErr)
  res |> mutate(gene = gene, draw = ds$label)
}

all_data <- lapply(draw_sizes, function(ds) {
  lapply(genes, load_results, ds = ds) |> bind_rows()
}) |> bind_rows()

ylims <- all_data |>
  mutate(lo = estimate - std_error, hi = estimate + std_error) |>
  group_by(draw) |>
  summarise(ymin = min(lo, na.rm = TRUE), ymax = max(hi, na.rm = TRUE),
            .groups = "drop") |>
  mutate(
    ymin = ymin - 0.005,
    ymax = ymax + 0.005
  )

# Force pct60 (row 4) to share the y range of pct30 (row 3).
ylims[ylims$draw == "pct60", c("ymin", "ymax")] <-
  ylims[ylims$draw == "pct30", c("ymin", "ymax")]


# ── Build panels ─────────────────────────────────────────────────────────────
# show_legend: show legend on the first panel of this column
# show_y:      show y-axis title, labels and ticks (FALSE for right column)
build_gene_column <- function(gene, show_legend = FALSE, show_y = TRUE) {
  panels <- list()

  for (i in seq_along(draw_sizes)) {
    ds      <- draw_sizes[[i]]
    dat     <- all_data |> filter(.data$gene == .env$gene, draw == ds$label)
    if (nrow(dat) == 0) {
      message("  [", gene, " ", ds$label, "] not found, skipping")
      next
    }

    ylim_row <- ylims |> filter(draw == ds$label)

    dat <- dat |>
      mutate(sig = classify_sig(p.value)) |>
      arrange(estimate) |>
      mutate(rank = row_number())

    is_first <- show_legend && i == 1

    p <- ggplot(dat, aes(x = rank, y = estimate, colour = sig)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.4) +
      geom_errorbar(
        aes(ymin = estimate - std_error, ymax = estimate + std_error),
        width = 0, alpha = 0.35, linewidth = 0.2
      ) +
      geom_point(size = 0.6, alpha = 0.75) +
      annotate("text", x = -Inf, y = Inf,
               label = paste0(gene, " — ", ds$label, " — ",
                              nrow(dat), " samples of ",
                              unique(dat$n), " sequences each"),
               hjust = -0.05, vjust = 1.4,
               size = 7 / .pt, colour = "grey20") +
      scale_colour_manual(
        name   = NULL,
        values = sig_colours,
        drop   = FALSE,
        guide  = guide_legend(
          override.aes = list(size = 2, alpha = 1),
          nrow = 1
        )
      ) +
      scale_y_continuous(limits = c(ylim_row$ymin, ylim_row$ymax)) +
      labs(
        x = "Sample rank (ordered by estimate)",
        y = "Temperature coefficient"
      ) +
      theme_bw(base_size = 9) +
      theme(
        panel.grid.minor     = element_blank(),
        panel.grid.major.x   = element_blank(),
        legend.position      = if (is_first) c(0.98, 0.05) else "none",
        legend.justification = c(1, 0),
        legend.direction     = "horizontal",
        legend.background    = element_rect(fill = "transparent", colour = NA),
        legend.text          = element_text(size = 7),
        legend.key.size      = unit(0.45, "lines")
      )

    panels[[ds$label]] <- p
  }

  if (length(panels) == 0) stop("No jackknife results found for ", gene)

  # Strip x-axis from all panels except the last (bottom row)
  last_label <- tail(names(panels), 1)
  panels <- lapply(names(panels), function(lbl) {
    p <- panels[[lbl]]
    if (lbl != last_label)
      p <- p + theme(axis.title.x = element_blank(),
                     axis.text.x  = element_blank(),
                     axis.ticks.x = element_blank())
    if (!show_y)
      p <- p + theme(axis.title.y = element_blank(),
                     axis.text.y  = element_blank(),
                     axis.ticks.y = element_blank())
    p
  })

  wrap_plots(panels, ncol = 1)
}

message("Building H1 column...")
col_h1 <- build_gene_column("H1", show_legend = TRUE,  show_y = TRUE)

message("Building H3 column...")
col_h3 <- build_gene_column("H3", show_legend = FALSE, show_y = FALSE)

fig <- col_h1 | col_h3

ggsave(file.path(out_dir, "supp3_jackknife.png"), fig,
       width = 200, height = 160, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp3_jackknife.pdf"), fig,
       width = 200, height = 160, units = "mm")

message("Written: supp3_jackknife.png / .pdf")
