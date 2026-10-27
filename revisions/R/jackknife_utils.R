# revisions/R/jackknife_utils.R
# Reading jackknife result tables.
#
# The revised jackknife scripts (03_rates_analysis/07_phylolm_jackknife.R and
# 04_epitope_analysis/21_jackknife_slac.R) fit the FULL model — six climate
# predictors plus the year factor — and store EVERY coefficient, one row per
# replicate per term. That means any predictor can be plotted later without
# refitting, but it also means the files are N_replicates x N_terms rows rather
# than N_replicates: 1000 replicates of a six-predictor model is 6000 rows.
#
# The submitted scripts fitted a temperature-only formula and stored one row per
# replicate, with no `term` column. read_jackknife() handles both, so a plotting
# script works against old and new files alike and a row count is never mistaken
# for a replicate count.

#' Read a jackknife results file, reduced to one coefficient.
#'
#' @param path Path to a *_jackknife_results.csv / *_results.csv file.
#' @param term Coefficient to keep. Ignored for legacy files, which contain only
#'   the temperature coefficient and carry no `term` column.
#' @param quiet Suppress the format message.
#' @return Data frame with one row per replicate.
read_jackknife <- function(path, term = "scaled_temp", quiet = FALSE) {
  if (!file.exists(path)) stop("Jackknife results not found: ", path)
  df <- readr::read_csv(path, show_col_types = FALSE)

  if (!"term" %in% names(df)) {
    if (!quiet)
      message("  ", basename(path), ": legacy format (", nrow(df),
              " rows, temperature only)")
    return(df)
  }

  avail <- unique(df$term)
  if (!term %in% avail)
    stop("Coefficient '", term, "' not in ", basename(path),
         ". Available: ", paste(avail, collapse = ", "))

  out <- df[df$term == term, , drop = FALSE]
  if (!quiet)
    message("  ", basename(path), ": ", nrow(df), " rows = ",
            nrow(out), " replicates x ", length(avail), " terms; kept '",
            term, "'")
  out
}

#' Summarise jackknife replicates for one coefficient.
#'
#' @param df Output of read_jackknife().
#' @return One-row data frame: replicate count, median, 95% range, and the
#'   share of replicates significant at p < 0.05 and sharing the median's sign.
summarise_jackknife <- function(df) {
  data.frame(
    replicates    = nrow(df),
    n_per_draw    = if ("n" %in% names(df)) df$n[1] else NA_integer_,
    median        = stats::median(df$estimate, na.rm = TRUE),
    q025          = unname(stats::quantile(df$estimate, 0.025, na.rm = TRUE)),
    q975          = unname(stats::quantile(df$estimate, 0.975, na.rm = TRUE)),
    pct_p05       = 100 * mean(df$p.value < 0.05, na.rm = TRUE),
    pct_same_sign = 100 * mean(sign(df$estimate) ==
                               sign(stats::median(df$estimate, na.rm = TRUE)),
                               na.rm = TRUE)
  )
}
