# scripts/05_final_figs/supp3_jackknife_sequences.R
#
# Supplementary Figure 3: Per-sequence jackknife effect size distributions
# for H1 and H3. Each panel shows sequences ordered by median temperature
# coefficient, with IQR and 5-95th percentile ribbons coloured by how many
# times each sequence was sampled (shared colour scale across both panels).
#
# INPUT:
#   data/rates/jackknife/{H1,H3}/jackknife_results.csv
#   data/rates/jackknife/{H1,H3}/jackknife_sequence_membership.csv
#
# OUTPUT: paper/figures/supp3_jackknife_sequences.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS <- 9

load_seq_summary <- function(gene) {
  jk_dir    <- file.path(PATHS$rates, "jackknife", gene)
  results   <- read_csv(file.path(jk_dir, "jackknife_results.csv"),
                        show_col_types = FALSE)
  if ("StdErr" %in% names(results) && !"std_error" %in% names(results))
    results <- rename(results, std_error = StdErr)
  membership <- read_csv(
    file.path(jk_dir, "jackknife_sequence_membership.csv"),
    show_col_types = FALSE
  )
  n_runs <- nrow(results)
  membership |>
    left_join(results |> select(sample_id, estimate), by = "sample_id") |>
    group_by(accession) |>
    summarise(
      median_est   = median(estimate, na.rm = TRUE),
      q25          = quantile(estimate, 0.25, na.rm = TRUE),
      q75          = quantile(estimate, 0.75, na.rm = TRUE),
      q05          = quantile(estimate, 0.05, na.rm = TRUE),
      q95          = quantile(estimate, 0.95, na.rm = TRUE),
      sample_freq  = n() / n_runs * 100,
      .groups      = "drop"
    ) |>
    arrange(median_est) |>
    mutate(rank = row_number(), gene = gene)
}

summaries <- lapply(c("H1", "H3"), load_seq_summary)
global_freq_range <- range(sapply(summaries, function(s) range(s$sample_freq)))

make_seq_panel <- function(seq_summary, show_x = FALSE) {
  gene <- unique(seq_summary$gene)
  panel_label <- ifelse(gene == "H1", "H1N1", "H3N2")

  ggplot(seq_summary, aes(x = rank, fill = sample_freq, colour = sample_freq)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25, colour = NA) +
    geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.45, colour = NA) +
    geom_line(aes(y = median_est), linewidth = 0.5) +
    scale_fill_viridis_c(
      name   = "Sampling\nfrequency (%)",
      option = "D",
      limits = global_freq_range
    ) +
    scale_colour_viridis_c(
      name   = "Sampling\nfrequency (%)",
      option = "D",
      limits = global_freq_range
    ) +
    annotate("text", x = -Inf, y = Inf, label = panel_label,
             hjust = -0.05, vjust = 1.4,
             size = BS / .pt, fontface = "bold", colour = "grey20") +
    labs(x = if (show_x) "Sequence rank (ordered by median estimate)" else NULL,
         y = NULL) +
    theme_bw(base_size = BS) +
    theme(
      panel.grid.minor = element_blank(),
      axis.text.x      = if (show_x) element_text() else element_blank(),
      axis.ticks.x     = if (show_x) element_line() else element_blank(),
      plot.margin      = margin(2, 4, 1, 4, "mm")
    )
}

p_H1 <- make_seq_panel(summaries[[1]], show_x = FALSE)
p_H3 <- make_seq_panel(summaries[[2]], show_x = TRUE)

fig <- (p_H1 / p_H3) +
  plot_layout(guides = "collect") &
  theme(legend.position = "right")

# Wrap with shared y-axis title
fig_labelled <- wrap_elements(full = fig) +
  labs(tag = "Regression coefficient of temperature") +
  theme(
    plot.tag          = element_text(size = BS, angle = 90, hjust = 0.5,
                                     vjust = 1),
    plot.tag.position = c(0.01, 0.5)
  )

ggsave(file.path(out_dir, "supp3_jackknife_sequences.png"), fig_labelled,
       width = 180, height = 120, units = "mm", dpi = 300)
ggsave(file.path(out_dir, "supp3_jackknife_sequences.pdf"), fig_labelled,
       width = 180, height = 120, units = "mm")

message("Written: supp3_jackknife_sequences.png / .pdf")