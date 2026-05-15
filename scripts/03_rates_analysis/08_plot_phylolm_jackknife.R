# scripts/03_rates_analysis/08_plot_phylolm_jackknife.R
#
# PURPOSE: Visualise jackknife subsampling results from 07_phylolm_jackknife.R.
#
# PLOT 1: Effect size (scaled_temp coefficient) per sample (1-100), coloured
#         by significance level (ns / * / ** / ***).
#
# PLOT 2: Per-sequence ribbon plot of effect sizes across all samples in which
#         each sequence appeared, ordered by median estimate, coloured by
#         how many times each sequence was sampled.
#
# INPUT:
#   data/rates/jackknife/{gene}/jackknife_results.csv
#   data/rates/jackknife/{gene}/jackknife_sequence_membership.csv
#
# OUTPUT:
#   figures/jackknife/{gene}_jackknife_samples.pdf
#   figures/jackknife/{gene}_jackknife_samples_ordered.pdf
#   figures/jackknife/{gene}_jackknife_sequences.png
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
})

source(here("config.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  GENES
}

fig_dir <- here("figures", "jackknife")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

SIG_LEVELS  <- c("ns", "*", "**", "***")
SIG_COLOURS <- setNames(viridis(4, end = 0.92), SIG_LEVELS)

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

for (gene in genes) {
  message("\n[", gene, "] Plotting jackknife results...")

  jk_dir   <- file.path(PATHS$rates, "jackknife", gene)
  res_path <- file.path(jk_dir, "jackknife_results.csv")
  mem_path <- file.path(jk_dir, "jackknife_sequence_membership.csv")

  if (!file.exists(res_path)) {
    warning(
      "Results not found: ", res_path,
      " — run 07_phylolm_jackknife.R first"
    )
    next
  }

  results    <- read_csv(res_path, show_col_types = FALSE)
  if ("StdErr" %in% names(results) && !"std_error" %in% names(results))
    results <- rename(results, std_error = StdErr)
  membership <- read_csv(mem_path, show_col_types = FALSE)

  # ── Plot 1a/1b: Effect size per sample — by order and sorted ──────────────
  results <- results |>
    mutate(sig = classify_sig(p.value))

  plot_samples <- function(dat, x_var, x_lab) {
    ggplot(dat, aes(x = {{ x_var }}, y = estimate, colour = sig)) +
      geom_hline(
        yintercept = 0, linetype = "dashed",
        colour = "grey50", linewidth = 0.4
      ) +
      geom_errorbar(
        aes(ymin = estimate - std_error, ymax = estimate + std_error),
        width = 0.4, alpha = 0.6, linewidth = 0.35
      ) +
      geom_point(size = 1.5) +
      scale_colour_manual(name = "p (raw)", values = SIG_COLOURS) +
      labs(
        title    = paste0(gene, " — jackknife temperature coefficients"),
        subtitle = paste0(
          nrow(dat), " samples of ", unique(dat$n), " sequences each"
        ),
        x = x_lab,
        y = "Estimate (scaled_temp)"
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())
  }

  p1a <- plot_samples(results, sample_id, "Sample (by run order)")

  results_sorted <- results |>
    arrange(estimate) |>
    mutate(rank = row_number())
  p1b <- plot_samples(results_sorted, rank, "Sample (ordered by estimate)")

  for (ext in c(".pdf", ".png")) {
    ggsave(file.path(fig_dir, paste0(gene, "_jackknife_samples", ext)),
           p1a, width = 9, height = 4, dpi = 200)
    ggsave(file.path(fig_dir, paste0(gene, "_jackknife_samples_ordered", ext)),
           p1b, width = 9, height = 4, dpi = 200)
  }
  message("  Written: ", gene, "_jackknife_samples(.pdf/.png) and _ordered")

  # ── Plot 2: Per-sequence effect size distribution ──────────────────────────
  # Pre-summarise to median + IQR + 5-95th pctile — far faster than per-seq
  # boxplots. Colour encodes how many times each sequence was sampled.
  seq_summary <- membership |>
    left_join(
      results |> select(sample_id, estimate),
      by = "sample_id"
    ) |>
    group_by(accession) |>
    summarise(
      median_est = median(estimate, na.rm = TRUE),
      q25        = quantile(estimate, 0.25, na.rm = TRUE),
      q75        = quantile(estimate, 0.75, na.rm = TRUE),
      q05        = quantile(estimate, 0.05, na.rm = TRUE),
      q95        = quantile(estimate, 0.95, na.rm = TRUE),
      n_sampled  = n(),
      .groups    = "drop"
    ) |>
    arrange(median_est) |>
    mutate(rank = row_number())

  n_seq <- nrow(seq_summary)

  p2 <- ggplot(
    seq_summary,
    aes(x = rank, fill = n_sampled, colour = n_sampled)
  ) +
    geom_hline(
      yintercept = 0, linetype = "dashed",
      colour = "grey50", linewidth = 0.4
    ) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25, colour = NA) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.45, colour = NA) +
    geom_line(aes(y = median_est), linewidth = 0.5) +
    scale_fill_viridis_c(name = "Times\nsampled", option = "D") +
    scale_colour_viridis_c(name = "Times\nsampled", option = "D") +
    labs(
      title    = paste0(gene, " — per-sequence jackknife effect sizes"),
      subtitle = paste0(
        "n = ", n_seq, " sequences ordered by median estimate; ",
        "ribbons = IQR (dark) and 5-95th percentile (light)"
      ),
      x = "Sequence rank (ordered by median estimate)",
      y = "Estimate (scaled_temp)"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  out2 <- file.path(fig_dir, paste0(gene, "_jackknife_sequences.png"))
  ggsave(out2, p2, width = 14, height = 5, dpi = 200)
  message("  Written: ", out2)
}

message("\nDone.")
