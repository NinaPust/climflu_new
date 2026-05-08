# scripts/04_epitope_analysis/06_compute_epitope_rates.R
#
# PURPOSE: Compute per-tip epitope and non-epitope dN/dS rates by reading
#          the fastme.bal trees from 05_map_to_tree.R, extracting terminal
#          branch lengths, and dividing by timetree branch length.
#          Joins all new rate columns onto the existing rates CSV.
#
# INPUT:
#   data/epitope/trees/{gene}_{type}_{dN|dS}_tree.nwk
#   data/trees/timed/{gene}_timetree.nwk
#   data/rates/{gene}_rates.csv
#
# OUTPUT:
#   data/rates/{gene}_rates_epitope.csv
#
# FLAGS:
#   --gene X         run a single gene (default: H1 H3)
#   --cutoff         apply log-rate floors (default: FALSE)
#   --omega_floor X  floor for log(omega)   (default: -8)
#   --dN_floor X     floor for log(dN_rate) (default: -20)
#   --dS_floor X     floor for log(dS_rate) (default: -20)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
})

source(here("config.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  c("H1", "H3")
}

apply_cutoff <- "--cutoff" %in% args
omega_floor  <- if ("--omega_floor" %in% args) {
  as.numeric(args[which(args == "--omega_floor") + 1])
} else -8
dN_floor     <- if ("--dN_floor" %in% args) {
  as.numeric(args[which(args == "--dN_floor") + 1])
} else -20
dS_floor     <- if ("--dS_floor" %in% args) {
  as.numeric(args[which(args == "--dS_floor") + 1])
} else -20

if (apply_cutoff)
  message("Applying log-rate floors: omega=", omega_floor,
          "  dN=", dN_floor, "  dS=", dS_floor)

epi_dir  <- here("data", "epitope")
tree_dir <- file.path(epi_dir, "trees")
types    <- c("epitope", "nonepitope", "ridge", "nonridge")

tip_branch_lengths <- function(tree) {
  tip_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  bl <- tree$edge.length[tip_idx]
  bl[bl < 0] <- 0
  setNames(bl, tree$tip.label)
}

for (gene in genes) {
  message("\n[", gene, "] Computing epitope/non-epitope rates...")

  timetree_path <- file.path(PATHS$trees_timed,
                             paste0(gene, "_timetree.nwk"))
  rates_path    <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))
  out_path      <- file.path(PATHS$rates,
                             paste0(gene, "_rates_epitope.csv"))

  if (!file.exists(timetree_path)) {
    warning("Timetree not found: ", timetree_path); next
  }
  if (!file.exists(rates_path)) {
    warning("Rates CSV not found: ", rates_path); next
  }

  timetree   <- read.tree(timetree_path)
  time_bl    <- tip_branch_lengths(timetree)
  base_rates <- read_csv(rates_path, show_col_types = FALSE)
  new_cols   <- tibble(accession = base_rates$accession)

  for (type in types) {
    dn_path <- file.path(tree_dir,
                         paste0(gene, "_", type, "_dN_tree.nwk"))
    ds_path <- file.path(tree_dir,
                         paste0(gene, "_", type, "_dS_tree.nwk"))

    if (!file.exists(dn_path) || !file.exists(ds_path)) {
      message("  [", type, "] Trees not found, filling NA")
      new_cols[[paste0(type, "_dN_rate")]]  <- NA_real_
      new_cols[[paste0(type, "_dS_rate")]]  <- NA_real_
      new_cols[[paste0(type, "_omega")]]    <- NA_real_
      next
    }

    dn_bl <- tip_branch_lengths(read.tree(dn_path))
    ds_bl <- tip_branch_lengths(read.tree(ds_path))

    common <- intersect(names(dn_bl),
                        intersect(names(ds_bl), names(time_bl)))
    t_bl   <- time_bl[common]
    valid  <- t_bl > 0

    all_acc      <- new_cols$accession
    dN_rate      <- setNames(rep(NA_real_, length(all_acc)), all_acc)
    dS_rate      <- setNames(rep(NA_real_, length(all_acc)), all_acc)
    branch_omega <- setNames(rep(NA_real_, length(all_acc)), all_acc)

    hits <- intersect(common[valid], all_acc)
    dN_rate[hits]      <- dn_bl[hits] / time_bl[hits]
    dS_rate[hits]      <- ds_bl[hits] / time_bl[hits]
    branch_omega[hits] <- dn_bl[hits] / pmax(ds_bl[hits], 1e-10)

    new_cols[[paste0(type, "_dN_rate")]]  <- dN_rate
    new_cols[[paste0(type, "_dS_rate")]]  <- dS_rate
    new_cols[[paste0(type, "_omega")]]    <- branch_omega

    message("  [", type, "] ", sum(!is.na(dN_rate)), " tips with rates")
  }

  # Log-transform; optionally apply floors to remove near-zero noise
  for (col in setdiff(names(new_cols), "accession")) {
    log_col  <- paste0("log_", col)
    vals     <- log(new_cols[[col]])
    vals[!is.finite(vals)] <- NA
    if (apply_cutoff) {
      floor <- if (grepl("_omega$",   col)) omega_floor else
               if (grepl("_dN_rate$", col)) dN_floor    else
               if (grepl("_dS_rate$", col)) dS_floor    else NA_real_
      if (!is.na(floor)) {
        n_before <- sum(!is.na(vals))
        vals[!is.na(vals) & vals < floor] <- NA
        n_after  <- sum(!is.na(vals))
        message("    floor log(", col, ") < ", floor, ": ",
                n_before - n_after, " removed, ", n_after, " kept")
      }
    }
    new_cols[[log_col]] <- vals
  }

  # Z-score scale
  for (col in grep("^log_", names(new_cols), value = TRUE)) {
    new_cols[[paste0("scaled_", col)]] <- scale(new_cols[[col]])[, 1]
  }

  df_out <- left_join(base_rates, new_cols, by = "accession")
  write_csv(df_out, out_path)
  message("  Written: ", out_path)
}

message("\nDone.")
