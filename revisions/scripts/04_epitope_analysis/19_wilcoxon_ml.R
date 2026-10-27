# revisions/scripts/04_epitope_analysis/19_wilcoxon_ml.R
#
# PURPOSE: Test whether recurrent non-synonymous mutations are associated with
#          climate, using the maximum-likelihood mutation calls from
#          18_ml_mutations.R. Revised counterpart to 11_wilcoxon_climate.R,
#          which ran on parsimony calls.
#
#          The statistical design is unchanged from script 11, deliberately: the
#          reviewers did not object to it, and keeping it identical means any
#          difference in the results is attributable to the ancestral
#          reconstruction rather than to the test.
#
# DESIGN:
#   For each mutation occurring on more than --min-count terminal branches,
#   split the tips into those carrying it and those not, and compare their
#   climate values with a Wilcoxon rank-sum test. Effect size is the
#   rank-biserial correlation r, which runs from -1 to +1:
#
#     r > 0   tips carrying the mutation sit at HIGHER values of that predictor
#     r < 0   they sit at LOWER values
#     r ~ 0   no shift
#
#   BH-FDR is applied within each gene x region, over all mutation x predictor
#   tests, as in script 11.
#
#   Climate predictors are the RAW values (not the scaled ones used by the
#   regression scripts), so medians print in degrees Celsius and percent.
#
# NOTE ON INDEPENDENCE: this test treats tips as independent, which they are
#   not — they share ancestry. A mutation inherited by a clade that happens to
#   circulate in a warm region will look climate-associated. Script 11 had the
#   same property. It is defensible here because mutations are called on
#   TERMINAL branches only, so each event arose independently on its own branch
#   rather than being inherited; but clustered sampling of related lineages can
#   still induce correlation, and that belongs in the discussion.
#
# INPUT:
#   revisions/data/epitope/{gene}_{region}_mutations_ml.csv   from script 18
#   revisions/data/rates/{gene}_rates_slac.csv                climate columns
#     (falls back to {gene}_rates_epitope.csv, then data/rates/{gene}_rates.csv)
#
# OUTPUT:
#   revisions/data/rates/{region}_mutation_climate_results_ml.csv
#   revisions/figures/epitope/{region}/mutation_climate_forest_ml.png
#   revisions/figures/epitope/{region}/mutation_climate_heatmap_ml.png
#
# FLAGS:
#   --region X      ridge | epitope (default: both)
#   --gene X        run a single gene (default: H1 H3)
#   --min-count N   a mutation qualifies if it occurs on MORE than N terminal
#                   branches (default 5, matching script 11)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(forcats)
  library(viridisLite)
})

source(here("revisions", "config.R"))

args      <- commandArgs(trailingOnly = TRUE)
regions   <- if ("--region" %in% args) args[which(args == "--region") + 1] else
             c("ridge", "epitope")
genes     <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
min_count <- if ("--min-count" %in% args) as.integer(args[which(args == "--min-count") + 1]) else 5L

predictors <- c(
  era5_temp_mean = "Temperature",
  humidity_rh    = "Relative humidity",
  delta_temp     = "ΔTemperature",
  delta_humidity = "ΔHumidity"
)

#' Rank-biserial correlation from a two-sample Wilcoxon statistic.
#'
#' r = 2U/(n1*n2) - 1, bounded [-1, 1]. Positive means the first group sits at
#' higher values.
rank_biserial <- function(x_pos, x_neg) {
  n1 <- length(x_pos); n2 <- length(x_neg)
  if (n1 == 0 || n2 == 0) return(NA_real_)
  unname(2 * wilcox.test(x_pos, x_neg, exact = FALSE)$statistic / (n1 * n2) - 1)
}

#' Locate a rates table carrying the raw climate columns.
rates_for <- function(gene) {
  cand <- c(file.path(PATHS$rates, paste0(gene, "_rates_slac.csv")),
            file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv")),
            file.path(here(), "data", "rates", paste0(gene, "_rates.csv")))
  cand <- cand[file.exists(cand)]
  for (p in cand) {
    nm <- names(suppressWarnings(read_csv(p, n_max = 1, show_col_types = FALSE)))
    if (all(names(predictors) %in% nm)) return(p)
  }
  NULL
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for (region in regions) {
  message("\n[", region, "] Wilcoxon climate analysis (ML mutation calls)...")

  out_dir <- here("revisions", "figures", "epitope", region)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  results <- lapply(genes, function(g) {
    mut_path <- file.path(PATHS$epitope,
                          paste0(g, "_", region, "_mutations_ml.csv"))
    if (!file.exists(mut_path)) {
      warning("Not found: ", basename(mut_path), " — run 18_ml_mutations.R first")
      return(NULL)
    }
    rates_path <- rates_for(g)
    if (is.null(rates_path)) {
      warning(g, ": no rates table with columns ",
              paste(names(predictors), collapse = ", "))
      return(NULL)
    }

    mut <- read_csv(mut_path, show_col_types = FALSE)
    rates <- suppressWarnings(read_csv(rates_path, show_col_types = FALSE,
                                       guess_max = 50000)) |>
      select(accession, all_of(names(predictors)))

    qualifying <- mut |> count(change) |> filter(n > min_count) |> pull(change)
    message("  [", g, "] ", nrow(mut), " events, ", n_distinct(mut$change),
            " distinct; qualifying (>", min_count, "): ", length(qualifying))
    if (length(qualifying) == 0) return(NULL)

    # One row per tip per qualifying mutation; tips with none get NA.
    tip_mut <- rates |>
      left_join(mut |> filter(change %in% qualifying) |> select(accession, change),
                by = "accession")

    lapply(qualifying, function(mut_name) {
      pos <- tip_mut |> filter(!is.na(change), change == mut_name)
      neg <- tip_mut |> filter(is.na(change) | change != mut_name)
      lapply(names(predictors), function(pred) {
        x_pos <- pos[[pred]][is.finite(pos[[pred]])]
        x_neg <- neg[[pred]][is.finite(neg[[pred]])]
        if (length(x_pos) < 3) return(NULL)
        wt <- wilcox.test(x_pos, x_neg, exact = FALSE)
        tibble(
          gene            = g,
          subtype         = ifelse(g %in% c("H1", "N1"), "H1N1", "H3N2"),
          region          = region,
          mutation        = mut_name,
          n_mutation      = length(x_pos),
          n_other         = length(x_neg),
          predictor       = pred,
          predictor_label = unname(predictors[[pred]]),
          median_mut      = median(x_pos),
          median_other    = median(x_neg),
          effect_r        = rank_biserial(x_pos, x_neg),
          p.value         = wt$p.value
        )
      }) |> bind_rows()
    }) |> bind_rows() |>
      mutate(p.fdr = p.adjust(p.value, method = "BH"))
  }) |> bind_rows()

  if (is.null(results) || nrow(results) == 0) {
    message("  No results, skipping plots."); next
  }

  res_path <- file.path(PATHS$rates,
                        paste0(region, "_mutation_climate_results_ml.csv"))
  write_csv(results, res_path)
  message("  Written: ", basename(res_path), "  (", nrow(results), " tests)")

  sig <- results |> filter(!is.na(p.fdr), p.fdr < 0.05)
  if (nrow(sig) > 0) {
    message("  Significant after BH-FDR:")
    print(as.data.frame(
      sig |> mutate(effect_r = round(effect_r, 3),
                    median_mut = round(median_mut, 2),
                    median_other = round(median_other, 2),
                    p.fdr = signif(p.fdr, 3)) |>
        select(subtype, mutation, predictor_label, n_mutation,
               median_mut, median_other, effect_r, p.fdr) |>
        arrange(p.fdr)), row.names = FALSE)
  } else {
    message("  No tests significant after FDR correction.")
  }

  pred_colours <- setNames(viridis(length(predictors), end = 0.9),
                           unname(predictors))
  turbo3 <- viridis(3, option = "turbo")

  # ── Forest plot ─────────────────────────────────────────────────────────────
  plot_df <- results |>
    mutate(sig        = !is.na(p.fdr) & p.fdr < 0.05,
           mutation   = fct_reorder(mutation, effect_r, .fun = median),
           point_fill = ifelse(sig, pred_colours[predictor_label], "white"))

  p1 <- ggplot(plot_df, aes(x = effect_r, y = mutation,
                            colour = predictor_label, group = predictor_label)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    geom_linerange(aes(xmin = 0, xmax = effect_r), linewidth = 0.6,
                   position = position_dodge(width = 0.6)) +
    geom_point(aes(fill = point_fill), shape = 21, size = 2.5, stroke = 0.7,
               position = position_dodge(width = 0.6)) +
    geom_text(data = plot_df |> distinct(subtype, mutation, n_mutation),
              aes(x = -Inf, y = mutation, label = paste0("n=", n_mutation)),
              inherit.aes = FALSE, hjust = -0.1, size = 2.5, colour = "grey40") +
    facet_wrap(~ subtype, scales = "free_y") +
    scale_colour_manual(values = pred_colours, name = NULL) +
    scale_fill_identity() +
    scale_x_continuous(limits = c(-1, 1),
                       breaks = c(-0.5, -0.3, -0.1, 0, 0.1, 0.3, 0.5)) +
    labs(title    = paste0("Climate associations of recurrent ", region,
                           " mutations (ML calls)"),
         subtitle = paste0("Rank-biserial correlation  |  ",
                           "Filled = BH-FDR p < 0.05  |  n = events"),
         x = "Effect size (rank-biserial r)", y = NULL) +
    theme_bw(base_size = 12) +
    theme(legend.position    = "bottom",
          panel.grid.minor   = element_blank(),
          panel.grid.major.y = element_blank())

  ggsave(file.path(out_dir, "mutation_climate_forest_ml.png"), p1,
         width = 11, height = 8, dpi = 180)
  message("  Saved: ", region, "/mutation_climate_forest_ml.png")

  # ── Heatmap ─────────────────────────────────────────────────────────────────
  p2 <- results |>
    mutate(sig = !is.na(p.fdr) & p.fdr < 0.05) |>
    ggplot(aes(x = predictor_label,
               y = fct_reorder(mutation, effect_r, .fun = median),
               fill = effect_r)) +
    geom_tile(colour = "white", linewidth = 0.4) +
    geom_tile(data = ~ filter(.x, sig), colour = "black", linewidth = 0.8,
              fill = NA) +
    geom_text(aes(label = round(effect_r, 2)), size = 2.8, colour = "white") +
    facet_wrap(~ subtype, scales = "free_y") +
    scale_fill_gradient2(low = turbo3[1], mid = "grey92", high = turbo3[3],
                         midpoint = 0, limits = c(-1, 1), name = "Effect r") +
    labs(title    = paste0("Climate effect sizes, recurrent ", region,
                           " mutations (ML calls)"),
         subtitle = "Black outline = BH-FDR p < 0.05", x = NULL, y = NULL) +
    theme_bw(base_size = 12) +
    theme(axis.text.x      = element_text(angle = 30, hjust = 1),
          panel.grid       = element_blank(),
          strip.background = element_blank(),
          strip.text       = element_text(face = "bold"))

  ggsave(file.path(out_dir, "mutation_climate_heatmap_ml.png"), p2,
         width = 10, height = 8, dpi = 180)
  message("  Saved: ", region, "/mutation_climate_heatmap_ml.png")
}

message("\nDone.")
