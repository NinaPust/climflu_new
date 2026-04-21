# R/climate_utils.R
# Utilities for joining climate data and scaling environmental predictors.

#' Join climate data onto a metadata data frame by accession.
#'
#' @param metadata  Data frame with an `accession` column.
#' @param climate   Data frame with an `accession` column and climate fields.
#' @return Joined data frame.
join_climate <- function(metadata, climate) {
  dplyr::left_join(metadata, climate, by = "accession")
}


#' Log-transform and z-score scale environmental predictors in place.
#'
#' Adds columns: log_rate, log_omega, scaled_temp, scaled_humidity,
#' and (if present) scaled_pop_density.
#'
#' @param df Data frame with rate, omega, era5_temp_mean, humidity_rh columns.
#' @return Data frame with additional scaled columns.
scale_env_vars <- function(df) {
  if ("rate" %in% names(df)) {
    df$log_rate <- log(df$rate)
    df$log_rate[!is.finite(df$log_rate)] <- NA
  }
  if ("omega" %in% names(df)) {
    df$log_omega <- log(df$omega)
    df$log_omega[!is.finite(df$log_omega)] <- NA
  }
  if ("era5_temp_mean" %in% names(df)) {
    df$scaled_temp <- scale(df$era5_temp_mean)[, 1]
  }
  if ("humidity_rh" %in% names(df)) {
    df$scaled_humidity <- scale(df$humidity_rh)[, 1]
  }
  if ("pop_density" %in% names(df)) {
    df$scaled_pop_density <- scale(log1p(df$pop_density))[, 1]
  }
  if (all(c("month", "latitude") %in% names(df))) {
    # Hemisphere-adjusted cosine of seasonal cycle.
    # Shift SH months by 6 so that peak flu season aligns with month 1 in both
    # hemispheres. Result: +1 = peak flu season, -1 = trough.
    adj_month        <- ifelse(df$latitude >= 0, df$month,
                               ((df$month - 1L + 6L) %% 12L) + 1L)
    df$season_cos    <- cos(2 * pi * (adj_month - 1) / 12)
    df$scaled_season <- scale(df$season_cos)[, 1]
  }
  df
}
