# revisions/R/forest_utils.R
# Shared forest-plot builder for the hurdle-model scripts (12a, 12b, 12c).
#
# Kept here rather than duplicated per script so the three analyses are visually
# comparable: same colour mapping per predictor, same filled/hollow convention
# for significance, same facet layout.

#' Standard predictor labels and viridis colours.
#'
#' Order matters — it fixes the dodge order and the legend order across all
#' three scripts.
PREDICTOR_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative humidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔHumidity",
  scaled_sampling_density = "Sampling density",
  scaled_pop_density      = "Population density"
)

PREDICTOR_COLOURS <- stats::setNames(
  viridisLite::viridis(length(PREDICTOR_LABELS), end = 0.92),
  PREDICTOR_LABELS
)

#' Region display order. Combined with `scale_y_discrete(limits = rev)` in the
#' plot functions, this puts the whole gene at the TOP and narrows down to the
#' ridge at the bottom, so a panel reads from least to most specific region.
REGION_LEVELS <- c("Full gene", "HA1", "HA2", "Epitope", "Ridge")

#' Region key -> display label. HA1 here means HA1 EXCLUDING the antigenic
#' sites, because the site map relabels those as epitope_*.
REGION_LABELS <- c(ridge = "Ridge", epitope = "Epitope", HA1 = "HA1",
                   HA2 = "HA2", full = "Full gene")


#' Prepare a results table for plotting.
#'
#' Adds 95% CIs, maps terms and regions to display labels, and derives subtype.
#'
#' @param res          Results data frame (from run_phyloglm / run_phyloglm_gee).
#' @param region_from  Column holding the region key (default "subset").
#' @return Data frame filtered to climate predictors, with plotting columns.
prep_forest <- function(res, region_from = "subset") {
  res <- res[res$term %in% names(PREDICTOR_LABELS), , drop = FALSE]
  if (nrow(res) == 0) stop("No climate predictor rows to plot")

  res$ci_lo     <- res$estimate - 1.96 * res$std_error
  res$ci_hi     <- res$estimate + 1.96 * res$std_error
  res$predictor <- factor(PREDICTOR_LABELS[res$term], levels = PREDICTOR_LABELS)
  res$region    <- factor(REGION_LABELS[res[[region_from]]],
                          levels = REGION_LEVELS)
  res$subtype   <- ifelse(res$gene %in% c("H1", "N1"), "H1N1", "H3N2")
  res[!is.na(res$region), , drop = FALSE]
}


#' Build a forest plot of coefficients by region and predictor.
#'
#' Filled points mark significance on `sig_col`; hollow points are n.s. Error
#' bars are 95% Wald intervals. Faceted by subtype so H1N1 and H3N2 sit side by
#' side on a shared x-scale.
#'
#' @param df       Output of prep_forest().
#' @param sig_col  Logical column name driving fill ("sig_raw" or "sig_fdr").
#' @param title    Plot title.
#' @param subtitle Plot subtitle.
#' @param xlab     x-axis label (coefficient scale differs by model).
#' @return A ggplot object.
plot_forest <- function(df, sig_col = "sig_fdr",
                        title = NULL, subtitle = NULL,
                        xlab = "Coefficient estimate") {
  df$point_fill <- ifelse(df[[sig_col]],
                          PREDICTOR_COLOURS[as.character(df$predictor)],
                          "white")

  ggplot2::ggplot(df, ggplot2::aes(x = estimate, y = region,
                                   colour = predictor, group = predictor)) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                        colour = "grey50", linewidth = 0.4) +
    ggplot2::geom_linerange(ggplot2::aes(xmin = ci_lo, xmax = ci_hi),
                            linewidth = 0.7,
                            position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::geom_point(ggplot2::aes(fill = point_fill), shape = 21,
                        size = 2.8, stroke = 0.8,
                        position = ggplot2::position_dodge(width = 0.7)) +
    ggplot2::scale_colour_manual(values = PREDICTOR_COLOURS, name = NULL) +
    ggplot2::scale_fill_identity() +
    ggplot2::scale_y_discrete(limits = rev, drop = TRUE) +
    ggplot2::facet_wrap(~ subtype) +
    ggplot2::labs(title = title, subtitle = subtitle, x = xlab, y = NULL) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      strip.text         = ggplot2::element_text(face = "bold", size = 11),
      strip.background   = ggplot2::element_blank(),
      legend.position    = "bottom",
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank()
    ) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2))
}


#' Write both the raw-p and FDR versions of a forest plot.
#'
#' Two figures rather than one because the raw-p version shows what the models
#' found and the FDR version shows what survives correction across all
#' gene x region x predictor tests — reviewers generally want to see both.
#'
#' @param res      Results data frame.
#' @param out_dir  Output directory.
#' @param stem     Filename stem; files are <stem>.png and <stem>_fdr.png.
#' @param title    Plot title.
#' @param xlab     x-axis label.
#' @param width,height,dpi ggsave dimensions.
#' @return Invisibly, the character vector of files written.
save_forest_pair <- function(res, out_dir, stem, title,
                             xlab = "Coefficient estimate",
                             width = 10, height = 6, dpi = 180) {
  df <- prep_forest(res)
  df$sig_raw <- !is.na(df$p.value) & df$p.value < 0.05
  df$sig_fdr <- !is.na(df$p.fdr)   & df$p.fdr   < 0.05

  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  files <- character(0)

  for (v in list(list(col = "sig_raw", sfx = "",
                      sub = "Coefficients ± 95% CI  |  Filled = p < 0.05"),
                 list(col = "sig_fdr", sfx = "_fdr",
                      sub = "Coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05"))) {
    p <- plot_forest(df, sig_col = v$col, title = title,
                     subtitle = v$sub, xlab = xlab)
    f <- file.path(out_dir, paste0(stem, v$sfx, ".png"))
    ggplot2::ggsave(f, p, width = width, height = height, dpi = dpi)
    message("  Saved: ", basename(f))
    files <- c(files, f)
  }
  invisible(files)
}
