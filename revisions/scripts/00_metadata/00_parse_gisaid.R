# scripts/00_metadata/00_parse_gisaid.R
#
# PURPOSE: Parse manually exported GISAID EpiFlu metadata XLS files into
#          standardised CSVs for downstream analysis. Also merges split FASTA
#          files (H3/N2/PA_H3N2 were downloaded in two batches: 2000-18 + 2019).
#
# NOTE: This script is subject to change as metadata processing methods improve.
#
# INPUT:
#   data/raw/H1N1_meta.xls                 — GISAID isolate metadata for H1N1
#   data/raw/H3N2_2000_18_meta.xls         — GISAID isolate metadata for H3N2 (2000-2018)
#   data/raw/H3N2_2019_metadata.xls        — GISAID isolate metadata for H3N2 (2019)
#   data/raw/{gene}_*.fasta                — GISAID sequence exports (merged here if split)
#
# OUTPUT:
#   data/raw/{gene}_sequences.fasta        — merged FASTA (canonical filename for pipeline)
#   data/metadata/{gene}_metadata.csv      — cleaned metadata with columns:
#     accession, strain, subtype, lineage, passage_history,
#     collection_date, year, month, location, country, region, host
#
# NOTE: lat/lon are not present in GISAID EpiFlu XLS exports. Coordinates will
#       be added in a later geocoding step.
#
# USAGE:
#   Rscript scripts/00_metadata/00_parse_gisaid.R
#   Rscript scripts/00_metadata/00_parse_gisaid.R --gene H3

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(readxl)
  library(stringr)
  library(lubridate)
})

source(here("revisions", "config.R"))

# ── Gene → raw file mapping ───────────────────────────────────────────────────
# FASTA parts to concatenate (in order). Output: {gene}_sequences.fasta
FASTA_PARTS <- list(
  H1      = c("H1_2000_19.fasta"),
  H3      = c("H3_2000_18.fasta", "H3_2019.fasta"),
  N1      = c("N1_2000_19.fasta"),
  N2      = c("N2_2000_18.fasta", "N2_2019.fasta"),
  PA_H1N1 = c("PA_H1N1_2000_19.fasta"),
  PA_H3N2 = c("PA_H3N2_2000_18.fasta", "PA_H3N2_2019.fasta")
)

# Metadata XLS files per gene (multiple = merge before processing).
# PA segments share isolate-level metadata with their parent subtype; the parse
# function filters metadata rows to accessions actually present in the FASTA.
META_FILES <- list(
  H1      = c("H1N1_meta.xls"),
  H3      = c("H3N2_2000_18_meta.xls", "H3N2_2019_metadata.xls"),
  N1      = c("H1N1_meta.xls"),
  N2      = c("H3N2_2000_18_meta.xls", "H3N2_2019_metadata.xls"),
  PA_H1N1 = c("H1N1_meta.xls"),
  PA_H3N2 = c("H3N2_2000_18_meta.xls", "H3N2_2019_metadata.xls")
)

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

# ── Parse command-line args ───────────────────────────────────────────────────
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else ALL_GENES

# ── Helper: merge FASTA parts ─────────────────────────────────────────────────
merge_fasta <- function(gene) {
  parts   <- FASTA_PARTS[[gene]]
  out_path <- file.path(PATHS$raw, paste0(gene, "_sequences.fasta"))

  if (length(parts) == 1 && parts[1] == paste0(gene, "_sequences.fasta")) {
    # Already the canonical name — nothing to do
    return(invisible(NULL))
  }

  part_paths <- file.path(PATHS$raw, parts)
  missing    <- part_paths[!file.exists(part_paths)]
  if (length(missing) > 0) {
    warning("Missing FASTA part(s), skipping merge for ", gene, ": ",
            paste(missing, collapse = ", "))
    return(invisible(NULL)) 
  }

  # Copy first part, then append any additional parts (OS-level, fast)
  file.copy(part_paths[1], out_path, overwrite = TRUE)
  if (length(part_paths) > 1) {
    for (p in part_paths[-1]) file.append(out_path, p)
  }

  n_seq <- sum(startsWith(read_lines(out_path, progress = FALSE), ">"))
  message("  FASTA merged: ", paste(parts, collapse = " + "),
          " → ", basename(out_path), " (", n_seq, " sequences)")
  invisible(out_path)
}

# ── Helper: parse one gene ────────────────────────────────────────────────────
parse_gisaid_metadata <- function(gene) {
  message("\n[", gene, "] Parsing GISAID metadata...")

  # Merge FASTA parts if needed
  merge_fasta(gene)

  # Read and combine metadata XLS file(s)
  xls_files <- file.path(PATHS$raw, META_FILES[[gene]])
  missing   <- xls_files[!file.exists(xls_files)]
  if (length(missing) > 0) {
    warning("Metadata file(s) not found, skipping ", gene, ": ",
            paste(missing, collapse = ", "))
    return(invisible(NULL))
  }

  raw <- bind_rows(lapply(xls_files, read_excel, col_types = "text"))

  # Deduplicate isolates that appear in both files (H3N2 overlap possible)
  raw <- distinct(raw, Isolate_Id, .keep_all = TRUE)

  # Filter metadata to accessions present in the FASTA. This is essential for
  # PA segments, whose metadata files cover all isolates of the subtype but
  # only a subset have PA sequences. It also guards against metadata/FASTA
  # mismatches for the other genes.
  fasta_path <- file.path(PATHS$raw, paste0(gene, "_sequences.fasta"))
  if (file.exists(fasta_path)) {
    all_lines <- read_lines(fasta_path, progress = FALSE)
    fasta_ids <- sub("^>", "", all_lines[startsWith(all_lines, ">")])
    n_fasta   <- length(fasta_ids)
    raw       <- filter(raw, Isolate_Id %in% fasta_ids)
    message("  Kept ", nrow(raw), " of ", n_fasta,
            " FASTA accessions found in metadata")
  } else {
    message("  Warning: FASTA not found at ", fasta_path,
            " — skipping accession filter")
  }

  message("  Rows after dedup + FASTA filter: ", nrow(raw))

  # Standardise to downstream column names
  df <- tibble(
    accession       = raw[["Isolate_Id"]],
    strain          = raw[["Isolate_Name"]],
    subtype         = raw[["Subtype"]],
    lineage         = raw[["Lineage"]],
    passage_history = raw[["Passage_History"]],
    collection_date = raw[["Collection_Date"]],
    location        = str_remove(raw[["Location"]], '^"+'),
    host            = raw[["Host"]]
  )

  # Parse collection date (GISAID uses YYYY-MM-DD, YYYY-MM, or YYYY)
  df <- df |>
    mutate(
      collection_date = parse_date_time(collection_date,
                                        orders = c("Y-m-d", "Y-m", "Y"),
                                        quiet = TRUE),
      year  = year(collection_date),
      month = month(collection_date)
    )

  # Parse country / region from location string ("Region / Country / ...")
  df <- df |>
    mutate(
      region  = str_split_fixed(location, " / ", 3)[, 1] |> str_trim(),
      country = str_split_fixed(location, " / ", 3)[, 2] |> str_trim()
    )

  # Coordinates are not present in EpiFlu XLS exports; set to NA for now
  df$latitude  <- NA_real_
  df$longitude <- NA_real_

  # Filter to date range and require accession
  n_before <- nrow(df)
  df <- df |>
    filter(
      !is.na(accession),
      !is.na(collection_date),
      collection_date >= as.Date(PARAMS$date_from),
      collection_date <= as.Date(PARAMS$date_to)
    )
  n_after <- nrow(df)
  message("  Rows: ", n_before, " → ", n_after,
          " after filtering (missing accession/date, out of range ",
          PARAMS$date_from, " – ", PARAMS$date_to, ")")

  # Reorder columns (lat/lon placeholders included for pipeline compatibility)
  df <- df |>
    select(accession, strain, subtype, lineage, passage_history,
           collection_date, year, month,
           location, country, region,
           latitude, longitude,
           host)

  out_path <- file.path(PATHS$metadata, paste0(gene, "_metadata.csv"))
  write_csv(df, out_path)
  message("  Written: ", out_path)

  invisible(df)
}

# ── Run ───────────────────────────────────────────────────────────────────────
results <- lapply(genes, parse_gisaid_metadata)
message("\nDone.")
