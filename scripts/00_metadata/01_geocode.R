# scripts/00_metadata/01_geocode.R
#
# PURPOSE: Geocode unique GISAID location strings using Nominatim (OpenStreetMap)
#          via the tidygeocoder package. Results are cached so re-runs only query
#          new locations. Coordinates are written back into all gene metadata CSVs.
#
# INPUT:   data/metadata/{gene}_metadata.csv  (run 00_parse_gisaid.R first)
# OUTPUT:  data/raw/location_coords_cache.csv — persistent geocoding cache
#          data/raw/geocode_failures.txt       — locations that could not be resolved
#          data/metadata/{gene}_metadata.csv   — lat/lon columns filled in
#
# RATE LIMIT: Nominatim requires <= 1 request/sec. tidygeocoder enforces this
#             automatically. For ~1000 unique locations expect ~20 min.
#
# USAGE:
#   Rscript scripts/00_metadata/01_geocode.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(tidygeocoder)
})

source(here("config.R"))

ALL_GENES    <- c(GENES, "PA_H1N1", "PA_H3N2")
cache_path   <- here("data", "raw", "location_coords_cache.csv")
failures_path <- here("data", "raw", "geocode_failures.txt")

max_fallback_levels <- 3   # how many location components to strip before giving up

# ── Step 1: collect unique location strings ───────────────────────────────────
message("Loading metadata and collecting unique locations...")

all_locations <- character(0)
for (g in ALL_GENES) {
  path <- file.path(PATHS$metadata, paste0(g, "_metadata.csv"))
  if (!file.exists(path)) {
    warning("Missing metadata, skipping: ", path)
    next
  }
  locs <- read_csv(path, col_select = "location",
                   show_col_types = FALSE)[["location"]]
  all_locations <- c(all_locations, locs)
}

unique_locs <- unique(all_locations[!is.na(all_locations) & all_locations != ""])
message("  Unique location strings: ", length(unique_locs))

# ── Step 2: load or initialise cache ─────────────────────────────────────────
if (file.exists(cache_path)) {
  cache <- read_csv(cache_path, show_col_types = FALSE)
  message("  Loaded cache: ", nrow(cache), " previously geocoded locations")
} else {
  cache <- tibble(
    location      = character(),
    lat           = numeric(),
    lon           = numeric(),
    geocode_level = integer(),   # 0 = exact match, 1+ = fallback levels used
    query_used    = character()
  )
  message("  No cache found — starting fresh")
}

to_geocode <- setdiff(unique_locs, cache$location)
message("  Locations to geocode: ", length(to_geocode))

# ── Step 3: geocode with fallback ─────────────────────────────────────────────

# Build a query string from a location string at a given fallback level.
# Drops the continent (level 1) and the N least-specific components from the end,
# then reverses the remaining parts so the most specific comes first.
# Example: "Africa / Kenya / Western Province / Kakamega County", level 0
#   → "Kakamega County, Western Province, Kenya"
# Example: same, level 1 (drop most specific sub-level)
#   → "Western Province, Kenya"
build_query <- function(location, drop_n = 0) {
  parts <- str_split(location, " / ")[[1]] |> str_trim()
  parts <- parts[-1]                          # drop continent
  if (drop_n > 0) parts <- head(parts, -drop_n)
  if (length(parts) == 0) return(NA_character_)
  paste(rev(parts), collapse = ", ")
}

geocode_one <- function(location) {
  for (level in 0:max_fallback_levels) {
    query <- build_query(location, drop_n = level)
    if (is.na(query)) break

    result <- geo(address = query, method = "osm",
                  lat = "lat", long = "lon",
                  quiet = TRUE)

    if (!is.na(result$lat)) {
      return(tibble(
        location      = location,
        lat           = result$lat,
        lon           = result$lon,
        geocode_level = as.integer(level),
        query_used    = query
      ))
    }
  }

  # All levels failed
  return(tibble(
    location      = location,
    lat           = NA_real_,
    lon           = NA_real_,
    geocode_level = NA_integer_,
    query_used    = NA_character_
  ))
}

if (length(to_geocode) > 0) {
  message("\nGeocoding ", length(to_geocode), " locations via Nominatim ",
          "(~1 req/sec, this may take a while)...")

  new_results <- vector("list", length(to_geocode))
  for (i in seq_along(to_geocode)) {
    new_results[[i]] <- geocode_one(to_geocode[i])
    if (i %% 50 == 0) {
      message("  ", i, " / ", length(to_geocode), " done")
    }
  }
  new_results <- bind_rows(new_results)

  n_ok   <- sum(!is.na(new_results$lat))
  n_fail <- sum(is.na(new_results$lat))
  message("  Resolved: ", n_ok, "  Failed: ", n_fail)

  # ── Step 4: update cache ──────────────────────────────────────────────────
  cache <- bind_rows(cache, new_results)
  write_csv(cache, cache_path)
  message("  Cache updated: ", nrow(cache), " total entries → ", cache_path)

  # Log failures
  failures <- new_results$location[is.na(new_results$lat)]
  if (length(failures) > 0) {
    writeLines(failures, failures_path)
    message("  Failures logged to: ", failures_path)
  }
} else {
  message("  All locations already cached — nothing to query.")
}

# ── Step 5: join coords back to all gene metadata CSVs ────────────────────────
message("\nWriting coordinates back to metadata CSVs...")

coords_lookup <- cache |>
  select(location, lat, lon) |>
  filter(!is.na(lat))

for (g in ALL_GENES) {
  path <- file.path(PATHS$metadata, paste0(g, "_metadata.csv"))
  if (!file.exists(path)) next

  df <- read_csv(path, show_col_types = FALSE)

  df <- df |>
    select(-latitude, -longitude) |>
    left_join(coords_lookup, by = "location") |>
    rename(latitude = lat, longitude = lon)

  # Restore column order expected downstream
  df <- df |>
    select(accession, strain, subtype, lineage, passage_history,
           collection_date, year, month,
           location, country, region,
           latitude, longitude,
           host)

  write_csv(df, path)
  n_coords <- sum(!is.na(df$latitude))
  message("  [", g, "] ", n_coords, " / ", nrow(df),
          " rows with coordinates → ", basename(path))
}

message("\nDone.")
