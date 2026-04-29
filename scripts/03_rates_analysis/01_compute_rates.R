# scripts/03_rates_analysis/01_compute_rates.R
#
# PURPOSE: Extract per-branch substitution rates from the time-calibrated tree
#          and join with dN/dS values and climate data.
#
# INPUT:
#   data/trees/timed/{gene}_timetree.nwk
#   data/trees/timed/{gene}_treetime_out/divergence_tree.nexus  (subst/site branch lengths)
#   data/trees/timed/{gene}_treetime_out/dates.tsv              (decimal dates per tip + root)
#   data/rates/{gene}_dNdS_matrix.rds
#   data/climate/{gene}_climate.csv
#   data/metadata/{gene}_metadata.csv
#
# OUTPUT:
#   data/rates/{gene}_rates.csv  — one row per tip, columns:
#     accession, rate (subst/site/yr), dN, dS, omega (dN/dS),
#     dN_rate, dS_rate, branch_omega  (per-branch rates from fastme trees),
#     era5_temp_mean, humidity_rh,
#     year, month, latitude, longitude, country, region, host,
#     log_rate, log_omega, log_dN_rate, log_dS_rate, log_branch_omega,
#     scaled_temp, scaled_humidity, scaled_pop_density
#
# FLAGS:
#   --gene X   run a single gene segment
#   --test     use subsampled data from data/test/

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
})

source(here("config.R"))
source(here("R", "climate_utils.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args     <- commandArgs(trailingOnly = TRUE)
genes    <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}
use_test <- "--test" %in% args

if (use_test) message("*** Running in TEST mode (subsampled data) ***")

paths_timed   <- if (use_test) PATHS$test_trees_timed else PATHS$trees_timed
paths_climate <- if (use_test) PATHS$test_climate     else PATHS$climate
paths_meta    <- if (use_test) PATHS$test_metadata    else PATHS$metadata
paths_rates   <- if (use_test) PATHS$test_rates       else PATHS$rates
paths_pop     <- PATHS$population

dir.create(paths_rates, recursive = TRUE, showWarnings = FALSE)

for (gene in genes) {
  message("\n[", gene, "] Computing rates...")

  timetree_path <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  dndS_path     <- file.path(PATHS$rates,   paste0(gene, "_dNdS_matrix.rds"))
  climate_path  <- file.path(paths_climate, paste0(gene, "_climate.csv"))
  meta_path     <- file.path(paths_meta,    paste0(gene, "_metadata.csv"))
  out_path      <- file.path(paths_rates,   paste0(gene, "_rates.csv"))

  missing <- c(timetree_path, climate_path, meta_path)[
    !file.exists(c(timetree_path, climate_path, meta_path))
  ]

  if (length(missing) > 0) {
    warning("Required file(s) not found, skipping ", gene, ":\n",
            paste(" ", missing, collapse = "\n"))
    next
  }

  treetime_dir  <- file.path(paths_timed, paste0(gene, "_treetime_out"))
  div_tree_path <- file.path(treetime_dir, "divergence_tree.nexus")

  if (!file.exists(div_tree_path))
    stop("Required TreeTime output not found: ", div_tree_path)

  divtree  <- read.nexus(div_tree_path)
  timetree <- read.tree(timetree_path)
  climate  <- read_csv(climate_path, show_col_types = FALSE)
  meta     <- read_csv(meta_path,    show_col_types = FALSE)

  has_dndS <- file.exists(dndS_path)
  if (has_dndS) {
    dndS <- readRDS(dndS_path)
  } else {
    message("  dNdS matrix not found — dN/dS columns will be NA")
  }

  # Terminal branch dN/dS rates from fastme.bal trees (GENES only, not PA)
  fastme_dir  <- file.path(PATHS$trees, "fastmebal")
  dn_tree_path <- file.path(fastme_dir, paste0(gene, "_dN_tree.nwk"))
  ds_tree_path <- file.path(fastme_dir, paste0(gene, "_dS_tree.nwk"))
  has_fastme <- file.exists(dn_tree_path) && file.exists(ds_tree_path)
  if (has_fastme) {
    dn_tree <- read.tree(dn_tree_path)
    ds_tree <- read.tree(ds_tree_path)
    tip_idx_dn <- match(seq_len(Ntip(dn_tree)), dn_tree$edge[, 2])
    tip_idx_ds <- match(seq_len(Ntip(ds_tree)), ds_tree$edge[, 2])
    dn_bl <- setNames(dn_tree$edge.length[tip_idx_dn], dn_tree$tip.label)
    ds_bl <- setNames(ds_tree$edge.length[tip_idx_ds], ds_tree$tip.label)
    dn_bl[dn_bl < 0] <- 0
    ds_bl[ds_bl < 0] <- 0
  } else {
    message("  fastme dN/dS trees not found — branch rate columns will be NA")
  }

  # ── Terminal branch rate (subst/site/yr) ──────────────────────────────────
  # divergence_tree.nexus: branch lengths in substitutions/site
  # timetree.nwk:          branch lengths in years (strict clock)
  # rate = terminal_branch_subst / terminal_branch_years
  # This gives the local evolutionary rate at the tip's sampling location,
  # not the cumulative root-to-tip average.
  div_tip_idx  <- match(seq_len(Ntip(divtree)),  divtree$edge[, 2])
  time_tip_idx <- match(seq_len(Ntip(timetree)), timetree$edge[, 2])

  div_bl  <- setNames(divtree$edge.length[div_tip_idx],   divtree$tip.label)
  time_bl <- setNames(timetree$edge.length[time_tip_idx], timetree$tip.label)

  common_tips <- intersect(names(div_bl), names(time_bl))
  common_tips <- common_tips[time_bl[common_tips] > 0]   # drop zero-time branches

  rate <- div_bl[common_tips] / time_bl[common_tips]
  message("  Clock rate from TreeTime: ",
          readLines(file.path(treetime_dir, "molecular_clock.txt"))[2])
  message("  Terminal branch rate — median: ",
          round(median(rate, na.rm = TRUE), 5),
          "  range: [", round(min(rate, na.rm = TRUE), 5),
          ", ", round(max(rate, na.rm = TRUE), 5), "]")

  # ── Mean dN, dS, omega per tip (row means of pairwise matrices) ───────────
  if (has_dndS) {
    common_acc <- intersect(names(rate), dndS$accessions)
    dN_mean    <- rowMeans(dndS$dN[common_acc, common_acc], na.rm = TRUE)
    dS_mean    <- rowMeans(dndS$dS[common_acc, common_acc], na.rm = TRUE)
    omega_mean <- dN_mean / dS_mean
  } else {
    common_acc <- names(rate)
    dN_mean    <- rep(NA_real_, length(common_acc))
    dS_mean    <- rep(NA_real_, length(common_acc))
    omega_mean <- rep(NA_real_, length(common_acc))
  }

  # ── Per-branch dN/dS rates (terminal branch dN or dS / time) ─────────────
  if (has_fastme) {
    branch_common <- intersect(common_acc, intersect(names(dn_bl), names(time_bl)))
    dN_rate    <- rep(NA_real_, length(common_acc))
    dS_rate    <- rep(NA_real_, length(common_acc))
    branch_omega <- rep(NA_real_, length(common_acc))
    names(dN_rate) <- names(dS_rate) <- names(branch_omega) <- common_acc
    t_bl <- time_bl[branch_common]
    valid <- t_bl > 0
    dN_rate[branch_common[valid]]    <- dn_bl[branch_common[valid]] / t_bl[valid]
    dS_rate[branch_common[valid]]    <- ds_bl[branch_common[valid]] / t_bl[valid]
    branch_omega[branch_common[valid]] <-
      dn_bl[branch_common[valid]] / pmax(ds_bl[branch_common[valid]], 1e-10)
  } else {
    dN_rate      <- rep(NA_real_, length(common_acc))
    dS_rate      <- rep(NA_real_, length(common_acc))
    branch_omega <- rep(NA_real_, length(common_acc))
  }

  # ── Assemble data frame ───────────────────────────────────────────────────
  df <- tibble(
    accession    = common_acc,
    rate         = rate[common_acc],
    dN           = dN_mean,
    dS           = dS_mean,
    omega        = omega_mean,
    dN_rate      = dN_rate,
    dS_rate      = dS_rate,
    branch_omega = branch_omega
  )

  pop_path <- file.path(paths_pop, paste0(gene, "_popdens.csv"))
  has_pop  <- file.exists(pop_path)
  if (has_pop) {
    pop <- read_csv(pop_path, show_col_types = FALSE) |>
      select(accession, pop_density, pop_epoch)
  } else {
    message("  Pop density file not found — pop_density will be NA")
  }

  df <- df |>
    left_join(climate |> select(accession, era5_temp_mean,
                                era5_dewpoint, humidity_rh,
                                any_of(c("delta_temp", "delta_humidity"))),
              by = "accession") |>
    left_join(meta |> select(accession, year, month, latitude, longitude,
                             country, region, host),
              by = "accession")

  if (has_pop) df <- left_join(df, pop, by = "accession")

  # ── Sampling density: sequences per (lat, lon, year, month) ──────────────
  samp_counts <- df |>
    filter(!is.na(latitude), !is.na(longitude), !is.na(year), !is.na(month)) |>
    count(latitude, longitude, year, month, name = "sampling_density")
  df <- left_join(df, samp_counts,
                  by = c("latitude", "longitude", "year", "month"))

  # ── Scale and transform predictors ───────────────────────────────────────
  df <- scale_env_vars(df)

  n_before <- nrow(df)
  df <- df |> filter(!is.na(rate), !is.na(era5_temp_mean), !is.na(humidity_rh))
  message("  Rows: ", n_before, " → ", nrow(df),
          " after dropping NA rates/climate")

  write_csv(df, out_path)
  message("  Written: ", out_path)
}

message("\nDone.")
