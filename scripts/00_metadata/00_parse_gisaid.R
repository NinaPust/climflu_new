# scripts/00_metadata/00_parse_gisaid.R
#
# PURPOSE: Parse manually exported GISAID metadata TSV files into standardised
#          CSVs for downstream analysis.
#
# NOTE: This script is subject to change as metadata processing methods improve.
#
# INPUT:
#   data/raw/{gene}_metadata.tsv  — GISAID metadata export (one per gene segment)
#   data/raw/{gene}_sequences.fasta — GISAID sequences (used only for accession check)
#
# OUTPUT:
#   data/metadata/{gene}_metadata.csv  — cleaned metadata with columns:
#     accession, strain, collection_date, year, month, latitude, longitude,
#     country, region, host
#
# GISAID EXPORT FIELDS REQUIRED (select when downloading):
#   Accession ID, Virus name, Collection date, Location,
#   Latitude, Longitude, Host
#
# USAGE:
#   Rscript scripts/00_metadata/00_parse_gisaid.R
#   Rscript scripts/00_metadata/00_parse_gisaid.R --gene H1

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(stringr)
  library(lubridate)
})

source(here("config.R"))
source(here("R", "climate_utils.R"))

# ── Parse command-line args ───────────────────────────────────────────────────
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

# ── Helper: parse one gene ────────────────────────────────────────────────────
parse_gisaid_metadata <- function(gene) {
  message("\n[", gene, "] Parsing GISAID metadata...")

  tsv_path <- file.path(PATHS$raw, paste0(gene, "_metadata.tsv"))
  if (!file.exists(tsv_path)) {
    warning("File not found, skipping: ", tsv_path)
    return(invisible(NULL))
  }

  raw <- read_tsv(tsv_path, show_col_types = FALSE)

  # Flexible column name mapping — GISAID exports vary slightly
  col_map <- list(
    accession       = c("Accession ID", "accession_id", "accession"),
    strain          = c("Virus name", "virus_name", "strain"),
    collection_date = c("Collection date", "collection_date", "date"),
    location        = c("Location", "location"),
    latitude        = c("Latitude", "latitude", "lat"),
    longitude       = c("Longitude", "longitude", "lon", "long"),
    host            = c("Host", "host")
  )

  resolve_col <- function(df, candidates) {
    match <- intersect(candidates, names(df))
    if (length(match) == 0) return(rep(NA_character_, nrow(df)))
    df[[match[1]]]
  }

  df <- tibble(
    accession       = resolve_col(raw, col_map$accession),
    strain          = resolve_col(raw, col_map$strain),
    collection_date = resolve_col(raw, col_map$collection_date),
    location        = resolve_col(raw, col_map$location),
    latitude        = as.numeric(resolve_col(raw, col_map$latitude)),
    longitude       = as.numeric(resolve_col(raw, col_map$longitude)),
    host            = resolve_col(raw, col_map$host)
  )

  # Parse date
  df <- df |>
    mutate(
      collection_date = parse_date_time(collection_date,
                                        orders = c("Y-m-d", "Y-m", "Y"),
                                        quiet = TRUE),
      year  = year(collection_date),
      month = month(collection_date)
    )

  # Parse country / region from location string (format: "Region / Country / ...")
  df <- df |>
    mutate(
      region  = str_split_fixed(location, " / ", 3)[, 1] |> str_trim(),
      country = str_split_fixed(location, " / ", 3)[, 2] |> str_trim()
    ) |>
    select(-location)

  # Filter
  n_before <- nrow(df)
  df <- df |>
    filter(
      !is.na(accession),
      !is.na(latitude), !is.na(longitude),
      !is.na(collection_date),
      collection_date >= as.Date(PARAMS$date_from),
      collection_date <= as.Date(PARAMS$date_to)
    )
  n_after <- nrow(df)
  message("  Rows: ", n_before, " → ", n_after, " after filtering (missing coords/dates, out-of-range)")

  out_path <- file.path(PATHS$metadata, paste0(gene, "_metadata.csv"))
  write_csv(df, out_path)
  message("  Written: ", out_path)

  invisible(df)
}

# ── Run ───────────────────────────────────────────────────────────────────────
results <- lapply(genes, parse_gisaid_metadata)
message("\nDone.")
