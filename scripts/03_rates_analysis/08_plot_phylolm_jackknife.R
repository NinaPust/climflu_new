# scripts/03_rates_analysis/08_plot_phylolm_jackknife.R
#
# PURPOSE: Visualise jackknife subsampling results from 07_phylolm_jackknife.R
#          for all four draw sizes (n1000, pct15, pct30, pct60).
#          Produces individual figures per draw size and one overview figure
#          with the samples-ordered plot stacked for all draw sizes.
#
# INPUT (per draw size):
#   data/rates/jackknife/{gene}/jackknife_results.csv           (n1000)
#   data/rates/jackknife/{gene}/{size}/jackknife_results.csv    (pct*)
#
# OUTPUT (per gene × draw size):
#   figures/jackknife/{gene}_{size}_jackknife_samples.pdf/.png
#   figures/jackknife/{gene}_{size}_jackknife_samples_ordered.pdf/.png
#   figures/jackknife/{gene}_{size}_jackknife_sequences.png
# OUTPUT (overview):
#   figures/jackknife/{gene}_jackknife_overview.pdf/.png
#
# FLAGS:
#   --gene X   run a single gene segment

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

args  <- commandArgs(trailingOnly = TRUE)
JK_COEF <- if ("--coef" %in% args) args[which(args == "--coef") + 1] else "scaled_temp"
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  GENES
}

draw_sizes <- list(
  list(label = "n1000", subdir = NULL),
  list(label = "pct15", subdir = "pct15"),
  list(label = "pct30", subdir = "pct30"),
  list(label = "pct60", subdir = "pct60")
)

sig_levels  <- c("ns", "*", "**", "***")
sig_colours <- setNames(viridis(4, end = 0.92), sig_levels)

classify_sig <- function(p) {
  factor(
    case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    levels = sig_levels
  )
}

for (gene in genes) {
  message("\n[", gene, "] Plotting jackknife overview...")

  fig_dir <- here("figures", "jackknife", gene)
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

  panels <- list()

  for (ds in draw_sizes) {
    jk_dir <- file.path(PATHS$rates, "jackknife", gene)
    if (!is.null(ds$subdir)) jk_dir <- file.path(jk_dir, ds$subdir)

    res_path <- file.path(jk_dir, "jackknife_results.csv")
    mem_path <- file.path(jk_dir, "jackknife_sequence_membership.csv")

    if (!file.exists(res_path)) {
      message("  [", ds$label, "] Results not found, skipping: ", res_path)
      next
    }

    message("  [", ds$label, "] Generating plots...")

    # Revised results store every coefficient (one row per replicate per
    # term), so reduce to the one being plotted. read_jackknife() also accepts
    # the submitted format, which had no `term` column.
    results <- read_jackknife(res_path, term = JK_COEF)
    if ("StdErr" %in% names(results) && !"std_error" %in% names(results)) {
      results <- rename(results, std_error = StdErr)
    }
    membership <- read_csv(mem_path, show_col_types = FALSE)

    results <- results |>
      mutate(sig = classify_sig(p.value))

    fig_prefix <- paste0(gene, "_", ds$label)

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
        scale_colour_manual(name = "p (raw)", values = sig_colours) +
        labs(
          title    = paste0(
            gene, " — jackknife temperature coefficients (", ds$label, ")"
          ),
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
      ggsave(
        file.path(fig_dir, paste0(fig_prefix, "_jackknife_samples", ext)),
        p1a, width = 9, height = 4, dpi = 200
      )
      ggsave(
        file.path(
          fig_dir, paste0(fig_prefix, "_jackknife_samples_ordered", ext)
        ),
        p1b, width = 9, height = 4, dpi = 200
      )
    }
    message(
      "  Written: ", fig_prefix, "_jackknife_samples(.pdf/.png) and _ordered"
    )

    # ── Plot 2: Per-sequence effect size distribution ────────────────────────
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
        title    = paste0(
          gene, " — per-sequence jackknife effect sizes (", ds$label, ")"
        ),
        subtitle = paste0(
          "n = ", n_seq, " sequences ordered by median estimate; ",
          "ribbons = IQR (dark) and 5-95th percentile (light)"
        ),
        x = "Sequence rank (ordered by median estimate)",
        y = "Estimate (scaled_temp)"
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())

    out2 <- file.path(fig_dir, paste0(fig_prefix, "_jackknife_sequences.png"))
    ggsave(out2, p2, width = 14, height = 5, dpi = 200)
    message("  Written: ", out2)

    # ── Collect ordered-samples panel for overview ───────────────────────────
    panels[[ds$label]] <- p1b +
      labs(
        title    = NULL,
        subtitle = paste0(
          ds$label, " — ",
          nrow(results_sorted), " samples of ",
          unique(results_sorted$n), " sequences each"
        )
      )
  }

  if (length(panels) == 0) {
    warning("  No results found for ", gene, ", skipping")
    next
  }

  overview <- wrap_plots(panels, ncol = 1) +
    plot_annotation(
      title = paste0(gene, " — jackknife temperature coefficients by draw size")
    )

  n_panels <- length(panels)
  for (ext in c(".pdf", ".png")) {
    ggsave(
      file.path(fig_dir, paste0(gene, "_jackknife_overview", ext)),
      overview, width = 9, height = 3.5 * n_panels, dpi = 200
    )
  }
  message(
    "  Written: ", gene, "_jackknife_overview (.pdf/.png) — ",
    n_panels, " panels"
  )
}

message("\nDone.")
