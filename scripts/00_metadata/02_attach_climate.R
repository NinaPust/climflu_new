# scripts/00_metadata/02_attach_climate.R
#
# PURPOSE: Fetch ERA5-Land monthly temperature and dewpoint from the
#          Copernicus Climate Data Store (CDS) for each sample's
#          lat/lon/date, compute relative humidity, and attach to
#          the metadata CSVs.
#
#          Downloads are done per unique year-month, extracting point
#          values for all locations in that month in one pass. Results
#          are cached so re-runs only fetch new months.
#
# SETUP (one-time):
#   1. Register at https://cds.climate.copernicus.eu
#   2. Copy your UID + API key from your CDS profile
#   3. Run in R:
#        ecmwfr::wf_set_key(user = "YOUR_UID", key = "YOUR_API_KEY",
#                            service = "cds")
#
# INPUT:
#   data/metadata/{gene}_metadata.csv  (with lat/lon from 01_geocode.R)
#
# OUTPUT:
#   data/climate/{gene}_climate.csv  — original metadata plus:
#       era5_temp_mean   (°C, monthly mean 2m temperature)
#       era5_temp_min    (°C, monthly minimum 2m temperature)
#       era5_temp_max    (°C, monthly maximum 2m temperature)
#       era5_dewpoint    (°C, monthly mean 2m dewpoint temperature)
#       humidity_rh      (%, relative humidity via Magnus formula)
#
# USAGE:
#   Rscript scripts/00_metadata/02_attach_climate.R
#   Rscript scripts/00_metadata/02_attach_climate.R --gene H1

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(terra)
  library(ecmwfr)
  library(ncdf4)
})

source(here("config.R"))

ALL_GENES   <- c(GENES, "PA_H1N1", "PA_H3N2")
cache_path  <- here("data", "climate", "era5_cache.csv")
nc_dir      <- here("data", "climate", "nc_downloads")
dir.create(nc_dir, recursive = TRUE, showWarnings = FALSE)

# ── Parse command-line args ───────────────────────────────────────────────────
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}

# ── Magnus formula for relative humidity ──────────────────────────────────────
relative_humidity <- function(temp_c, dewpoint_c) {
  a <- 17.625
  b <- 243.04
  gamma_t  <- (a * temp_c) / (b + temp_c)
  gamma_dp <- (a * dewpoint_c) / (b + dewpoint_c)
  rh <- 100.0 * exp(gamma_dp - gamma_t)
  round(rh, 2)
}

# ── Step 1: load metadata and find unique (lat, lon, year, month) ─────────────
message("Loading metadata...")

gene_dfs <- list()
for (g in genes) {
  path <- file.path(PATHS$metadata, paste0(g, "_metadata.csv"))
  if (!file.exists(path)) {
    warning("Missing metadata, skipping: ", path)
    next
  }
  gene_dfs[[g]] <- read_csv(path, show_col_types = FALSE)
}

if (length(gene_dfs) == 0) {
  stop("No metadata files found. Run 00_parse_gisaid.R and ",
       "01_geocode.R first.")
}

all_meta <- bind_rows(gene_dfs) |>
  filter(!is.na(latitude), !is.na(longitude),
         !is.na(year), !is.na(month))

unique_keys <- all_meta |>
  distinct(latitude, longitude, year, month) |>
  rename(lat = latitude, lon = longitude)

message("  Unique (lat, lon, year, month) combinations: ", nrow(unique_keys))

# ── Step 2: load or initialise cache ──────────────────────────────────────────
if (file.exists(cache_path)) {
  cache <- read_csv(cache_path, show_col_types = FALSE)
  message("  Loaded cache: ", nrow(cache), " entries")
} else {
  cache <- tibble(
    lat = numeric(), lon = numeric(),
    year = integer(), month = integer(),
    era5_temp_mean = numeric(), era5_temp_min = numeric(),
    era5_temp_max = numeric(), era5_dewpoint = numeric(),
    humidity_rh = numeric()
  )
}

# Filter to uncached keys
if (nrow(cache) > 0) {
  to_fetch <- unique_keys |>
    anti_join(cache, by = c("lat", "lon", "year", "month"))
} else {
  to_fetch <- unique_keys
}

message("  Already cached: ", nrow(unique_keys) - nrow(to_fetch))
message("  To fetch: ", nrow(to_fetch))

# ── Step 3: download ERA5-Land monthly data and extract point values ──────────

# Group uncached keys by year-month (one CDS request per month)
fetch_month <- function(yr, mn, points_df) {
  # points_df: tibble with lat, lon columns for this year-month
  nc_file <- file.path(nc_dir,
                       sprintf("era5land_%04d_%02d.nc", yr, mn))

  # Download if not already on disk
  if (!file.exists(nc_file)) {
    # Bounding box with 1° padding
    lat_min <- floor(min(points_df$lat)) - 1
    lat_max <- ceiling(max(points_df$lat)) + 1
    lon_min <- floor(min(points_df$lon)) - 1
    lon_max <- ceiling(max(points_df$lon)) + 1

    request <- list(
      product_type = "monthly_averaged_reanalysis",
      variable = c("2m_temperature", "2m_dewpoint_temperature"),
      year    = sprintf("%04d", yr),
      month   = sprintf("%02d", mn),
      time    = "00:00",
      area    = c(lat_max, lon_min, lat_min, lon_max),
      format  = "netcdf",
      dataset_short_name = "reanalysis-era5-land-monthly-means",
      target  = basename(nc_file)
    )

    message("    Downloading ERA5-Land for ", yr, "-",
            sprintf("%02d", mn), "...")

    wf_request(
      request = request,
      transfer = TRUE,
      path = nc_dir,
      verbose = FALSE
    )
  }

  # Extract point values from NetCDF
  r <- rast(nc_file)
  pts <- vect(points_df, geom = c("lon", "lat"), crs = "EPSG:4326")

  # Identify layers — ERA5-Land monthly means NetCDF typically has
  # layers named with the variable
  layer_names <- names(r)
  temp_layers <- grep("t2m|temperature_2m|2m_temperature",
                      layer_names, value = TRUE, ignore.case = TRUE)
  dew_layers  <- grep("d2m|dewpoint|2m_dewpoint",
                      layer_names, value = TRUE, ignore.case = TRUE)

  if (length(temp_layers) == 0 || length(dew_layers) == 0) {
    warning("Could not identify layers in ", nc_file,
            ". Layers: ", paste(layer_names, collapse = ", "))
    points_df$era5_temp_mean <- NA_real_
    points_df$era5_temp_min  <- NA_real_
    points_df$era5_temp_max  <- NA_real_
    points_df$era5_dewpoint  <- NA_real_
    points_df$humidity_rh    <- NA_real_
    return(points_df)
  }

  # Monthly averaged product gives one value per variable
  temp_vals <- terra::extract(r[[temp_layers[1]]], pts)[, 2]
  dew_vals  <- terra::extract(r[[dew_layers[1]]], pts)[, 2]

  # ERA5 stores temperature in Kelvin — convert to Celsius
  temp_c <- temp_vals - 273.15
  dew_c  <- dew_vals - 273.15

  # Monthly means product has one layer, so min/max = mean
  # (daily min/max require the daily product; for now use monthly mean)
  points_df$era5_temp_mean <- round(temp_c, 3)
  points_df$era5_temp_min  <- round(temp_c, 3)
  points_df$era5_temp_max  <- round(temp_c, 3)
  points_df$era5_dewpoint  <- round(dew_c, 3)
  points_df$humidity_rh    <- relative_humidity(temp_c, dew_c)

  points_df
}

if (nrow(to_fetch) > 0) {
  month_groups <- to_fetch |>
    mutate(ym = paste(year, month, sep = "-")) |>
    group_by(year, month) |>
    group_split()

  message("\nFetching ERA5-Land data for ",
          length(month_groups), " year-months...")

  new_results <- vector("list", length(month_groups))
  for (i in seq_along(month_groups)) {
    grp <- month_groups[[i]]
    yr  <- grp$year[1]
    mn  <- grp$month[1]
    message("  [", i, "/", length(month_groups), "] ",
            yr, "-", sprintf("%02d", mn),
            " (", nrow(grp), " points)")

    new_results[[i]] <- tryCatch(
      fetch_month(yr, mn, grp),
      error = function(e) {
        warning("Failed for ", yr, "-", mn, ": ", e$message)
        grp$era5_temp_mean <- NA_real_
        grp$era5_temp_min  <- NA_real_
        grp$era5_temp_max  <- NA_real_
        grp$era5_dewpoint  <- NA_real_
        grp$humidity_rh    <- NA_real_
        grp
      }
    )

    # Checkpoint cache every 10 months
    if (i %% 10 == 0) {
      cache <- bind_rows(cache, bind_rows(new_results[!sapply(new_results, is.null)]))
      write_csv(cache, cache_path)
      message("    Cache checkpoint: ", nrow(cache), " entries")
    }
  }

  new_results <- bind_rows(new_results)
  cache <- bind_rows(cache, new_results) |>
    distinct(lat, lon, year, month, .keep_all = TRUE)
  write_csv(cache, cache_path)
  message("  Cache updated: ", nrow(cache), " total entries")
}

# ── Step 4: join climate data back to metadata CSVs ───────────────────────────
message("\nWriting climate CSVs...")

coords_lookup <- cache |>
  select(lat, lon, year, month,
         era5_temp_mean, era5_temp_min, era5_temp_max,
         era5_dewpoint, humidity_rh)

for (g in names(gene_dfs)) {
  df <- gene_dfs[[g]]
  out_path <- file.path(PATHS$climate, paste0(g, "_climate.csv"))

  out_df <- df |>
    left_join(coords_lookup,
              by = c("latitude" = "lat", "longitude" = "lon",
                     "year" = "year", "month" = "month"))

  write_csv(out_df, out_path)
  n_ok <- sum(!is.na(out_df$era5_temp_mean))
  message("  [", g, "] ", n_ok, " / ", nrow(out_df),
          " rows with climate data → ", basename(out_path))
}

message("\nDone.")
