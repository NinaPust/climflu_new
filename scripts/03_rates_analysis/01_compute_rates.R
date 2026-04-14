# scripts/03_rates_analysis/01_compute_rates.R
#
# PURPOSE: Extract per-branch substitution rates from the time-calibrated tree
#          and join with dN/dS values and climate data.
#
# INPUT:
#   data/trees/{gene}_timetree.nwk
#   data/rates/{gene}_dNdS_matrix.rds
#   data/climate/{gene}_climate.csv
#   data/metadata/{gene}_metadata.csv
#
# OUTPUT:
#   data/rates/{gene}_rates.csv  — one row per tip, columns:
#     accession, rate (subst/site/yr), dN, dS, omega (dN/dS),
#     era5_temp_mean, era5_temp_min, humidity_rh,
#     year, month, latitude, longitude, country, region, host,
#     log_rate, log_omega, scaled_temp, scaled_humidity

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
})

source(here("config.R"))
source(here("R", "climate_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] Computing rates...")

  timetree_path <- file.path(PATHS$trees, paste0(gene, "_timetree.nwk"))
  dndS_path     <- file.path(PATHS$rates, paste0(gene, "_dNdS_matrix.rds"))
  climate_path  <- file.path(PATHS$climate, paste0(gene, "_climate.csv"))
  meta_path     <- file.path(PATHS$metadata, paste0(gene, "_metadata.csv"))
  out_path      <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))

  for (p in c(timetree_path, dndS_path, climate_path, meta_path)) {
    if (!file.exists(p)) stop("Required file not found: ", p)
  }

  timetree <- read.tree(timetree_path)
  dndS     <- readRDS(dndS_path)
  climate  <- read_csv(climate_path,  show_col_types = FALSE)
  meta     <- read_csv(meta_path,     show_col_types = FALSE)

  # ── Extract tip rates from timetree branch lengths ────────────────────────
  # Branch lengths in a time-calibrated tree are in years; dividing by the
  # root-to-tip distance gives substitution rate (subst/site/yr).
  # We use root-to-tip distance as a proxy for rate.
  rtt <- dist.nodes(timetree)
  root_idx <- Ntip(timetree) + 1L
  tip_idx  <- seq_len(Ntip(timetree))
  root_to_tip <- rtt[root_idx, tip_idx]
  names(root_to_tip) <- timetree$tip.label

  # Tree total time span (root to most recent tip)
  meta_dates <- meta |> filter(accession %in% timetree$tip.label)
  max_year   <- max(meta_dates$year, na.rm = TRUE)
  min_year   <- min(meta_dates$year, na.rm = TRUE)
  time_span  <- max_year - min_year
  if (time_span <= 0) stop("Cannot compute rates: tree spans <= 0 years")

  rate <- root_to_tip / time_span

  # ── Mean dN, dS, omega per tip (row means of pairwise matrices) ───────────
  common_acc <- intersect(names(rate), dndS$accessions)
  dN_mean  <- rowMeans(dndS$dN[common_acc, common_acc], na.rm = TRUE)
  dS_mean  <- rowMeans(dndS$dS[common_acc, common_acc], na.rm = TRUE)
  omega_mean <- dN_mean / dS_mean

  # ── Assemble data frame ───────────────────────────────────────────────────
  df <- tibble(
    accession = common_acc,
    rate      = rate[common_acc],
    dN        = dN_mean,
    dS        = dS_mean,
    omega     = omega_mean
  )

  df <- df |>
    left_join(climate, by = "accession") |>
    left_join(meta |> select(accession, year, month, latitude, longitude,
                              country, region, host),
              by = "accession")

  # ── Scale and transform predictors ───────────────────────────────────────
  df <- scale_env_vars(df)

  n_before <- nrow(df)
  df <- df |> filter(!is.na(rate), !is.na(era5_temp_mean), !is.na(humidity_rh))
  message("  Rows: ", n_before, " → ", nrow(df), " after dropping NA rates/climate")

  write_csv(df, out_path)
  message("  Written: ", out_path)
}

message("\nDone.")
