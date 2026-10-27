# scripts/02_phylogenetics/05_assemble_matrices.R
#
# PURPOSE: Stitch pairwise dN/dS matrices from yn00 into full distance
#          matrices and save assembled dN, dS, omega matrices.
#
# INPUT:
#   data/yn00/{gene}/  — directories of yn00 output
#   data/trees/pruned/{gene}_accessions.txt
#
# OUTPUT:
#   data/rates/{gene}_dNdS_matrix.rds   — full assembled matrices (dN, dS, omega)

suppressPackageStartupMessages({
  library(here)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "matrix_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] Assembling dN/dS matrices...")

  yn00_root <- file.path(PATHS$yn00, gene)
  acc_file  <- file.path(PATHS$metadata, paste0(gene, "_accessions_final.txt"))

  if (!dir.exists(yn00_root)) {
    warning("yn00 directory not found, skipping: ", yn00_root)
    next
  }

  # Accession list: prefer post-pruning list from trees/pruned, fallback to metadata
  pruned_acc <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
  acc_file   <- if (file.exists(pruned_acc)) pruned_acc else acc_file
  accessions <- readLines(acc_file)
  accessions <- accessions[nchar(trimws(accessions)) > 0]
  n <- length(accessions)
  message("  ", n, " accessions (from ", basename(acc_file), ")")

  # ── Collect pairwise matrices ─────────────────────────────────────────────
  pair_dirs <- list.dirs(yn00_root, full.names = TRUE, recursive = FALSE)
  pair_dirs <- pair_dirs[file.exists(file.path(pair_dirs, "2YN.dS"))]
  message("  ", length(pair_dirs), " completed codeml pairs found")

  # Prefer remapped RDS (full accession names); fall back to raw PAML matrices
  load_dn <- function(d) {
    rds <- file.path(d, "dN_remapped.rds")
    if (file.exists(rds)) readRDS(rds) else read_paml_matrix(file.path(d, "2YN.dN"))
  }
  load_ds <- function(d) {
    rds <- file.path(d, "dS_remapped.rds")
    if (file.exists(rds)) readRDS(rds) else read_paml_matrix(file.path(d, "2YN.dS"))
  }

  dN_matrices <- lapply(pair_dirs, load_dn)
  dS_matrices <- lapply(pair_dirs, load_ds)

  # ── Stitch into full matrices ─────────────────────────────────────────────
  dN_full <- matrix_stitch(dN_matrices, expected_labels = accessions)
  dS_full <- matrix_stitch(dS_matrices, expected_labels = accessions)

  # ── Validate and remove mismatched accessions ─────────────────────────────
  common <- intersect(rownames(dN_full), rownames(dS_full))
  common <- intersect(common, accessions)
  removed <- setdiff(accessions, common)
  if (length(removed) > 0) {
    message("  Removed ", length(removed), " accessions not present in both matrices:")
    message("    ", paste(head(removed, 10), collapse = ", "),
            if (length(removed) > 10) "..." else "")
  }

  dN_full <- dN_full[common, common]
  dS_full <- dS_full[common, common]

  validate_matrix(dN_full, expected_n = length(common))
  validate_matrix(dS_full, expected_n = length(common))

  # ── Save outputs ──────────────────────────────────────────────────────────
  saveRDS(
    list(dN = dN_full, dS = dS_full, omega = dN_full / dS_full, accessions = common),
    file.path(PATHS$rates, paste0(gene, "_dNdS_matrix.rds"))
  )

  message("  Matrices: ", file.path(PATHS$rates, paste0(gene, "_dNdS_matrix.rds")))
}

message("\nDone.")
