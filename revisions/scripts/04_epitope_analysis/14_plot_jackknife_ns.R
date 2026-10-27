# scripts/04_epitope_analysis/14_plot_jackknife_ns.R
#
# PURPOSE: Visualise jackknife results from 13_jackknife_ns_rates.R.
#          Mirrors scripts/03_rates_analysis/08_plot_phylolm_jackknife.R
#          but loops over regions and adds a combined summary panel.
#
# PLOT A: Effect size per sample (by run order), coloured by significance.
# PLOT B: Effect size per sample (ordered by estimate magnitude).
# PLOT C: Per-sequence ribbon plot ordered by median estimate.
# PLOT D: Combined faceted summary across all regions and genes.
#
# INPUT:
#   data/rates/jackknife/epitope/{gene}_{region}_jackknife_results.csv
#   data/rates/jackknife/epitope/{gene}_{region}_jackknife_membership.csv
#
# OUTPUT:
#   figures/jackknife/epitope/{gene}_{region}_jackknife_samples.pdf/png
#   figures/jackknife/epitope/{gene}_{region}_jackknife_samples_ordered.pdf/png
#   figures/jackknife/epitope/{gene}_{region}_jackknife_sequences.png
#   figures/jackknife/epitope/jackknife_ns_summary.png
#
# FLAGS:
#   --gene X     run a single gene (default: H1 H3)
#   --region X   run a single region (default: full HA2 epitope ridge)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
})

source(here("revisions", "config.R"))

args    <- commandArgs(trailingOnly = TRUE)
genes   <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
regions <- if ("--region" %in% args) args[which(args == "--region") + 1] else c("full", "HA2", "epitope", "ridge")

jk_dir  <- file.path(PATHS$rates, "jackknife", "epitope")
fig_dir <- here("revisions", "figures", "jackknife", "epitope")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

SIG_LEVELS  <- c("ns", "*", "**", "***")
SIG_COLOURS <- setNames(viridis(4, end = 0.92), SIG_LEVELS)

REGION_LABELS <- c(
  full    = "Whole gene",
  HA2     = "HA2",
  epitope = "Epitope",
  ridge   = "Ridge"
)

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

all_results <- list()

for (gene in genes) {
  for (region in regions) {
    tag      <- paste0("[", gene, " ", region, "]")
    res_path <- file.path(jk_dir, paste0(gene, "_", region,
                                          "_jackknife_results.csv"))
    mem_path <- file.path(jk_dir, paste0(gene, "_", region,
                                          "_jackknife_membership.csv"))

    if (!file.exists(res_path)) {
      warning(tag, " Results not found: ", res_path,
              " — run 13_jackknife_ns_rates.R first")
      next
    }

    message("\n", tag, " Plotting jackknife results...")

    results    <- read_csv(res_path, show_col_types = FALSE)
    membership <- read_csv(mem_path, show_col_types = FALSE)

    if ("StdErr" %in% names(results) && !"std_error" %in% names(results))
      results <- rename(results, std_error = StdErr)

    results <- results |> mutate(sig = classify_sig(p.value))
    all_results[[paste(gene, region, sep = "_")]] <- results

    stem  <- file.path(fig_dir, paste0(gene, "_", region, "_jackknife"))
    title <- paste0(gene, " ", REGION_LABELS[[region]],
                    " — jackknife temperature coefficients")
    sub   <- paste0(nrow(results), " samples of ",
                    unique(results$n), " sequences each")

    # ── Plot A: by run order ────────────────────────────────────────────────
    p_a <- ggplot(results, aes(x = sample_id, y = estimate, colour = sig)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.4) +
      geom_errorbar(aes(ymin = estimate - std_error,
                        ymax = estimate + std_error),
                    width = 0.4, alpha = 0.6, linewidth = 0.35) +
      geom_point(size = 1.5) +
      scale_colour_manual(name = "p (raw)", values = SIG_COLOURS) +
      labs(title = title, subtitle = sub,
           x = "Sample (by run order)", y = "Estimate (scaled_temp)") +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())

    for (ext in c(".pdf", ".png"))
      ggsave(paste0(stem, "_samples", ext), p_a, width = 9, height = 4, dpi = 200)

    # ── Plot B: ordered by estimate ─────────────────────────────────────────
    results_sorted <- results |>
      arrange(estimate) |>
      mutate(rank = row_number())

    p_b <- ggplot(results_sorted, aes(x = rank, y = estimate, colour = sig)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.4) +
      geom_errorbar(aes(ymin = estimate - std_error,
                        ymax = estimate + std_error),
                    width = 0.4, alpha = 0.6, linewidth = 0.35) +
      geom_point(size = 1.5) +
      scale_colour_manual(name = "p (raw)", values = SIG_COLOURS) +
      labs(title = title, subtitle = sub,
           x = "Sample (ordered by estimate)", y = "Estimate (scaled_temp)") +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())

    for (ext in c(".pdf", ".png"))
      ggsave(paste0(stem, "_samples_ordered", ext), p_b,
             width = 9, height = 4, dpi = 200)

    message("  Written: ", basename(stem), "_samples(.pdf/.png) and _ordered")

    # ── Plot C: per-sequence ribbon ─────────────────────────────────────────
    seq_summary <- membership |>
      left_join(results |> select(sample_id, estimate), by = "sample_id") |>
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

    p_c <- ggplot(seq_summary,
                  aes(x = rank, fill = n_sampled, colour = n_sampled)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey50", linewidth = 0.4) +
      geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25, colour = NA) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.45, colour = NA) +
      geom_line(aes(y = median_est), linewidth = 0.5) +
      scale_fill_viridis_c(name = "Times\nsampled", option = "D") +
      scale_colour_viridis_c(name = "Times\nsampled", option = "D") +
      labs(
        title    = paste0(gene, " ", REGION_LABELS[[region]],
                          " — per-sequence jackknife effect sizes"),
        subtitle = paste0(
          "n = ", nrow(seq_summary),
          " sequences ordered by median estimate; ",
          "ribbons = IQR (dark) and 5–95th percentile (light)"
        ),
        x = "Sequence rank (ordered by median estimate)",
        y = "Estimate (scaled_temp)"
      ) +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())

    ggsave(paste0(stem, "_sequences.png"), p_c,
           width = 14, height = 5, dpi = 200)
    message("  Written: ", basename(stem), "_sequences.png")
  }
}

# ── Combined summary panel ──────────────────────────────────────────────────
if (length(all_results) > 0) {
  message("\nBuilding combined summary panel...")

  summary_df <- bind_rows(all_results) |>
    mutate(
      subtype = ifelse(gene == "H1", "H1N1", "H3N2"),
      region  = factor(REGION_LABELS[region],
                       levels = REGION_LABELS[c("ridge", "epitope", "full", "HA2")])
    ) |>
    arrange(region, estimate) |>
    group_by(subtype, region) |>
    mutate(rank = row_number()) |>
    ungroup() |>
    mutate(
      ci_lo = estimate - 1.96 * std_error,
      ci_hi = estimate + 1.96 * std_error
    )

  # Per-facet sample-size annotation (only shown when n < 1000)
  annot_df <- summary_df |>
    group_by(subtype, region) |>
    summarise(
      draw_n   = unique(n)[1],
      pool_n   = round(unique(n)[1] / 0.6),
      .groups  = "drop"
    ) |>
    filter(draw_n < 1000) |>
    mutate(label = paste0("n per draw = ", draw_n,
                          "\n(60% of ", pool_n, ")"))

  p_summary <- ggplot(summary_df,
                      aes(x = rank, y = estimate, colour = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(ymin = ci_lo, ymax = ci_hi),
                   linewidth = 0.3, alpha = 0.4) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text(data = annot_df,
              aes(label = label),
              x = -Inf, y = Inf,
              hjust = -0.07, vjust = 1.3,
              inherit.aes = FALSE,
              size = 2.8, colour = "grey30") +
    scale_colour_manual(name = "p (raw)", values = SIG_COLOURS) +
    facet_grid(subtype ~ region) +
    labs(
      title = "Jackknife robustness: temperature effect on NS substitution rates",
      x     = "Sample (ordered by estimate)",
      y     = "Estimate (scaled_temp)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      panel.grid.minor  = element_blank(),
      strip.text        = element_text(face = "bold"),
      legend.position   = "bottom"
    )

  out_summary <- file.path(fig_dir, "jackknife_ns_summary.png")
  ggsave(out_summary, p_summary, width = 14, height = 7, dpi = 200)
  message("Written: jackknife_ns_summary.png")
}

message("\nDone.")
