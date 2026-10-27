# scripts/00_metadata/00z_explore_climate.R
#
# PURPOSE: Exploratory diagnostics on climate-attached metadata.
#          Gitignored — not part of the reproducible pipeline.
#
# INPUT:  data/climate/{gene}_climate.csv
# OUTPUT: figures/metadata_diagnostics/{gene}_climate_diagnostics.pdf
#
# USAGE:  Rscript scripts/00_metadata/00z_explore_climate.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(sf)
  library(rnaturalearth)
})

source(here("revisions", "config.R"))

out_dir <- here("revisions", "figures", "metadata_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes_plot <- c("H1", "H3", "N1", "N2", "PA_H1N1", "PA_H3N2")
world <- ne_countries(scale = "medium", returnclass = "sf")

for (g in genes_plot) {
  path <- file.path(PATHS$climate, paste0(g, "_climate.csv"))
  if (!file.exists(path)) {
    warning("Missing: ", path)
    next
  }

  df <- read_csv(path, show_col_types = FALSE) |>
    filter(!is.na(era5_temp_mean), !is.na(humidity_rh))

  if (nrow(df) == 0) {
    warning("No climate data for ", g)
    next
  }

  message("[", g, "] ", nrow(df), " rows with climate data")

  pdf_path <- file.path(out_dir, paste0(g, "_climate_diagnostics.pdf"))
  pdf(pdf_path, width = 10, height = 7)

  # ── 1. Temperature & humidity distributions ───────────────────────────
  p1a <- ggplot(df, aes(x = era5_temp_mean)) +
    geom_histogram(bins = 50, fill = "steelblue", colour = "white") +
    labs(title = paste(g, "— Monthly mean temperature"),
         x = "Temperature (°C)", y = "Count") +
    theme_minimal(base_size = 11)

  p1b <- ggplot(df, aes(x = humidity_rh)) +
    geom_histogram(bins = 50, fill = "darkorange", colour = "white") +
    labs(title = paste(g, "— Relative humidity"),
         x = "RH (%)", y = "Count") +
    theme_minimal(base_size = 11)

  print(p1a)
  print(p1b)

  # ── 2. Climate by top countries ───────────────────────────────────────
  top_countries <- df |>
    filter(!is.na(country), country != "") |>
    count(country, sort = TRUE) |>
    head(15) |>
    pull(country)

  p2 <- df |>
    filter(country %in% top_countries) |>
    mutate(country = factor(country, levels = rev(top_countries))) |>
    ggplot(aes(x = era5_temp_mean, y = country)) +
    geom_boxplot(fill = "steelblue", alpha = 0.6, outlier.size = 0.3) +
    labs(title = paste(g, "— Temperature by top 15 countries"),
         x = "Temperature (°C)", y = NULL) +
    theme_minimal(base_size = 11)

  print(p2)

  # ── 3. Seasonal climate by hemisphere ─────────────────────────────────
  df_season <- df |>
    filter(!is.na(month), !is.na(latitude)) |>
    mutate(
      hemisphere = ifelse(latitude >= 0, "Northern", "Southern"),
      month_label = factor(month, levels = 1:12, labels = month.abb)
    )

  p3a <- ggplot(df_season,
                aes(x = month_label, y = era5_temp_mean,
                    fill = hemisphere)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    scale_fill_manual(values = c("Northern" = "tomato",
                                 "Southern" = "dodgerblue")) +
    labs(title = paste(g, "— Temperature by month & hemisphere"),
         x = "Month", y = "Temperature (°C)", fill = NULL) +
    theme_minimal(base_size = 11)

  p3b <- ggplot(df_season,
                aes(x = month_label, y = humidity_rh,
                    fill = hemisphere)) +
    geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
    scale_fill_manual(values = c("Northern" = "tomato",
                                 "Southern" = "dodgerblue")) +
    labs(title = paste(g, "— Humidity by month & hemisphere"),
         x = "Month", y = "RH (%)", fill = NULL) +
    theme_minimal(base_size = 11)

  print(p3a)
  print(p3b)

  # ── 4. Map coloured by temperature ────────────────────────────────────
  p4 <- ggplot() +
    geom_sf(data = world, fill = "grey90", colour = "grey70",
            linewidth = 0.1) +
    geom_point(data = df |> filter(!is.na(latitude), !is.na(longitude)),
               aes(x = longitude, y = latitude, colour = era5_temp_mean),
               size = 0.4, alpha = 0.5) +
    scale_colour_viridis_c(option = "inferno") +
    labs(title = paste(g, "— Samples coloured by temperature"),
         colour = "Temp (°C)") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom", panel.grid = element_blank())

  print(p4)

  # ── 5. Temperature vs latitude ────────────────────────────────────────
  p5 <- df |>
    filter(!is.na(latitude)) |>
    ggplot(aes(x = latitude, y = era5_temp_mean)) +
    geom_point(size = 0.3, alpha = 0.2, colour = "steelblue") +
    geom_smooth(method = "lm", colour = "red", se = FALSE) +
    labs(title = paste(g, "— Temperature vs latitude (ERA5 sanity check)"),
         x = "Latitude", y = "Temperature (°C)") +
    theme_minimal(base_size = 11)

  print(p5)

  # ── 6. Temperature vs humidity ────────────────────────────────────────
  p6 <- ggplot(df, aes(x = era5_temp_mean, y = humidity_rh)) +
    geom_point(size = 0.3, alpha = 0.15, colour = "grey30") +
    geom_density_2d(colour = "tomato", linewidth = 0.4, bins = 8) +
    labs(title = paste(g, "— Temperature vs relative humidity"),
         x = "Temperature (°C)", y = "RH (%)") +
    theme_minimal(base_size = 11)

  print(p6)

  # ── 7–11. Delta plots (only if previous-month data available) ─────────
  has_delta <- all(c("delta_temp", "delta_humidity") %in% names(df)) &&
    any(!is.na(df$delta_temp))

  if (has_delta) {
    df_d <- df |> filter(!is.na(delta_temp), !is.na(delta_humidity))
    message("  ", nrow(df_d), " rows with delta climate data")

    # 7. Distributions of deltas
    p7a <- ggplot(df_d, aes(x = delta_temp)) +
      geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      labs(title = paste(g, "— Monthly temperature change (current − previous)"),
           x = "ΔTemp (°C)", y = "Count") +
      theme_minimal(base_size = 11)

    p7b <- ggplot(df_d, aes(x = delta_humidity)) +
      geom_histogram(bins = 60, fill = "darkorange", colour = "white") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "red") +
      labs(title = paste(g, "— Monthly humidity change (current − previous)"),
           x = "ΔRH (%)", y = "Count") +
      theme_minimal(base_size = 11)

    print(p7a)
    print(p7b)

    # 8. Seasonal pattern of deltas by hemisphere
    df_d_season <- df_d |>
      filter(!is.na(month), !is.na(latitude)) |>
      mutate(
        hemisphere  = ifelse(latitude >= 0, "Northern", "Southern"),
        month_label = factor(month, levels = 1:12, labels = month.abb)
      )

    p8a <- ggplot(df_d_season,
                  aes(x = month_label, y = delta_temp, fill = hemisphere)) +
      geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_fill_manual(values = c("Northern" = "tomato",
                                   "Southern" = "dodgerblue")) +
      labs(title = paste(g, "— ΔTemp by month & hemisphere"),
           x = "Month", y = "ΔTemp (°C)", fill = NULL) +
      theme_minimal(base_size = 11)

    p8b <- ggplot(df_d_season,
                  aes(x = month_label, y = delta_humidity, fill = hemisphere)) +
      geom_boxplot(outlier.size = 0.3, alpha = 0.7) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      scale_fill_manual(values = c("Northern" = "tomato",
                                   "Southern" = "dodgerblue")) +
      labs(title = paste(g, "— ΔHumidity by month & hemisphere"),
           x = "Month", y = "ΔRH (%)", fill = NULL) +
      theme_minimal(base_size = 11)

    print(p8a)
    print(p8b)

    # 9. Delta temp vs delta humidity (joint distribution)
    p9 <- ggplot(df_d, aes(x = delta_temp, y = delta_humidity)) +
      geom_point(size = 0.3, alpha = 0.15, colour = "grey30") +
      geom_density_2d(colour = "tomato", linewidth = 0.4, bins = 8) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50") +
      labs(title = paste(g, "— ΔTemp vs ΔHumidity"),
           x = "ΔTemp (°C)", y = "ΔRH (%)") +
      theme_minimal(base_size = 11)

    print(p9)

    # 10. Absolute temp vs delta temp (are large changes tied to level?)
    p10 <- ggplot(df_d, aes(x = era5_temp_mean, y = delta_temp)) +
      geom_point(size = 0.3, alpha = 0.15, colour = "steelblue") +
      geom_smooth(method = "lm", colour = "red", se = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      labs(title = paste(g, "— Temperature level vs monthly change"),
           x = "Temp (°C)", y = "ΔTemp (°C)") +
      theme_minimal(base_size = 11)

    p11 <- ggplot(df_d, aes(x = humidity_rh, y = delta_humidity)) +
      geom_point(size = 0.3, alpha = 0.15, colour = "darkorange") +
      geom_smooth(method = "lm", colour = "red", se = FALSE) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50") +
      labs(title = paste(g, "— Humidity level vs monthly change"),
           x = "RH (%)", y = "ΔRH (%)") +
      theme_minimal(base_size = 11)

    print(p10)
    print(p11)

  } else {
    message("  No delta columns found — skipping delta plots")
  }

  dev.off()
  message("  Saved: ", pdf_path)
}

message("\nAll climate diagnostics written to: ", out_dir)
