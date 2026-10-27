# scripts/04_epitope_analysis/11_wilcoxon_climate.R
#
# PURPOSE: Test whether specific recurrent NS mutations are associated with
#          climate conditions (Wilcoxon rank-sum + rank-biserial effect size).
#          Runs for any region that has a mutations CSV (ridge / epitope).
#
#          Qualifying mutations: > min_count independent occurrences.
#          BH-FDR correction within each gene × region.
#
#          Outputs per region to figures/epitope/{region}/:
#            mutation_climate_forest.png
#            mutation_climate_heatmap.png
#
# INPUT:
#   data/epitope/{gene}_{region}_mutations.csv   (from script 09)
#   data/rates/{gene}_rates_epitope.csv
#
# OUTPUT:
#   data/rates/{region}_mutation_climate_results.csv
#   figures/epitope/{region}/mutation_climate_forest.png
#   figures/epitope/{region}/mutation_climate_heatmap.png
#
# FLAGS:
#   --region X   ridge | epitope  (default: ridge epitope)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(viridisLite)
})

source(here("config.R"))

args    <- commandArgs(trailingOnly = TRUE)
regions <- if ("--region" %in% args) args[which(args == "--region") + 1] else
           c("ridge", "epitope")

min_count  <- 5L
predictors <- c(
  era5_temp_mean = "Temperature",
  humidity_rh    = "Relative humidity",
  delta_temp     = "ΔTemperature",
  delta_humidity = "ΔHumidity"
)

rank_biserial <- function(x_pos, x_neg) {
  n1 <- length(x_pos); n2 <- length(x_neg)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  2 * wilcox.test(x_pos, x_neg)$statistic / (n1 * n2) - 1
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for (region in regions) {
  message("\n[", region, "] Wilcoxon climate analysis...")

  out_dir <- here("figures", "epitope", region)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- lapply(c("H1", "H3"), function(g) {
    mut_path <- file.path(here("data", "epitope"),
                          paste0(g, "_", region, "_mutations.csv"))
    if (!file.exists(mut_path)) {
      warning("Mutations file not found: ", mut_path,
              "\nRun script 09 first.")
      return(NULL)
    }
    mut <- read_csv(mut_path, show_col_types = FALSE)
    rates <- read_csv(
      file.path(PATHS$rates, paste0(g, "_rates_epitope.csv")),
      show_col_types = FALSE
    ) |> select(accession, all_of(names(predictors)))

    qualifying <- mut |> count(change) |>
      filter(n > min_count) |> pull(change)

    message("  [", g, "] Qualifying mutations (>", min_count, "): ",
            length(qualifying))
    if (length(qualifying) == 0) return(NULL)

    tip_mut <- rates |>
      left_join(
        mut |> filter(change %in% qualifying) |> select(accession, change),
        by = "accession"
      )

    rows <- lapply(qualifying, function(mut_name) {
      pos <- tip_mut |> filter(change == mut_name)
      neg <- tip_mut |> filter(is.na(change) | change != mut_name)
      lapply(names(predictors), function(pred) {
        x_pos <- pos[[pred]][is.finite(pos[[pred]])]
        x_neg <- neg[[pred]][is.finite(neg[[pred]])]
        if (length(x_pos) < 3) return(NULL)
        wt <- wilcox.test(x_pos, x_neg, exact = FALSE)
        tibble(
          gene            = g,
          subtype         = ifelse(g == "H1", "H1N1", "H3N2"),
          mutation        = mut_name,
          n_mutation      = length(x_pos),
          n_other         = length(x_neg),
          predictor       = pred,
          predictor_label = predictors[[pred]],
          median_mut      = median(x_pos),
          median_other    = median(x_neg),
          effect_r        = rank_biserial(x_pos, x_neg),
          p.value         = wt$p.value
        )
      }) |> bind_rows()
    }) |> bind_rows()

    rows |> mutate(p.fdr = p.adjust(p.value, method = "BH"))
  }) |> bind_rows()

  if (is.null(results) || nrow(results) == 0) {
    message("  No results, skipping plots.")
    next
  }

  res_path <- file.path(PATHS$rates,
                        paste0(region, "_mutation_climate_results.csv"))
  write_csv(results, res_path)
  message("  Results written: ", res_path)

  sig <- results |> filter(!is.na(p.fdr), p.fdr < 0.05)
  if (nrow(sig) > 0) {
    message("  Significant after FDR (p.fdr < 0.05):")
    print(sig |> select(subtype, mutation, predictor_label, effect_r,
                        p.value, p.fdr) |> arrange(p.fdr))
  } else {
    message("  No significant results after FDR correction.")
  }

  pred_colours <- setNames(viridis(length(predictors), end = 0.9), predictors)
  turbo3       <- viridis(3, option = "turbo")

  # ── Forest plot ─────────────────────────────────────────────────────────────
  plot_df <- results |>
    mutate(
      sig        = !is.na(p.fdr) & p.fdr < 0.05,
      mutation   = fct_reorder(mutation, effect_r, .fun = median),
      point_fill = ifelse(sig, pred_colours[predictor_label], "white")
    )

  p1 <- ggplot(plot_df,
               aes(x = effect_r, y = mutation,
                   colour = predictor_label, group = predictor_label)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(xmin = 0, xmax = effect_r), linewidth = 0.6,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(fill = point_fill), shape = 21, size = 2.5, stroke = 0.7,
               position = position_dodge(width = 0.6)) +
    geom_text(
      data = plot_df |> distinct(subtype, mutation, n_mutation),
      aes(x = -Inf, y = mutation, label = paste0("n=", n_mutation)),
      inherit.aes = FALSE, hjust = -0.1, size = 2.5, colour = "grey40"
    ) +
    facet_wrap(~subtype, scales = "free_y") +
    scale_colour_manual(values = pred_colours, name = NULL) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(-1, 1),
                       breaks = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)) +
    labs(
      title    = paste0("Climate associations of recurrent ", region, " mutations"),
      subtitle = paste0("Rank-biserial correlation  |  ",
                        "Filled = BH-FDR p < 0.05  |  n = events"),
      x = "Effect size (rank-biserial r)", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position    = "bottom",
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank())

  ggsave(file.path(out_dir, "mutation_climate_forest.png"),
         p1, width = 11, height = 8, dpi = 180)
  message("  Saved: mutation_climate_forest.png")

  # ── Heatmap ──────────────────────────────────────────────────────────────────
  p2 <- results |>
    mutate(sig = !is.na(p.fdr) & p.fdr < 0.05) |>
    ggplot(aes(x = predictor_label,
               y = fct_reorder(mutation, effect_r, .fun = median),
               fill = effect_r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(data = ~ filter(.x, sig),
              colour = "black", linewidth = 0.8, fill = NA) +
    geom_text(aes(label = round(effect_r, 2)), size = 2.8, colour = "white") +
    facet_wrap(~subtype, scales = "free_y") +
    scale_fill_gradient2(low = turbo3[1], mid = "grey92", high = turbo3[3],
                         midpoint = 0, limits = c(-1, 1), name = "Effect r") +
    labs(
      title    = paste0("Climate effect sizes for recurrent ", region, " mutations"),
      subtitle = "Black outline = BH-FDR p < 0.05",
      x = NULL, y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(axis.text.x     = element_text(angle = 30, hjust = 1),
          legend.position = "right",
          panel.grid      = element_blank())

  ggsave(file.path(out_dir, "mutation_climate_heatmap.png"),
         p2, width = 11, height = 7, dpi = 180)
  message("  Saved: mutation_climate_heatmap.png")
}

message("\nDone.")
