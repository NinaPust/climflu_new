# scripts/02_phylogenetics/03_assemble_matrices.R
#
# PURPOSE: Stitch pairwise dN/dS matrices from codeml into full distance
#          matrices, validate dimensions, and build dN and dS trees.
#
# INPUT:
#   data/codeml/{gene}/  — directories of codeml output
#   data/metadata/{gene}_accessions_final.txt
#
# OUTPUT:
#   data/trees/{gene}_dN_tree.nwk
#   data/trees/{gene}_dS_tree.nwk
#   data/rates/{gene}_dNdS_matrix.rds   — full assembled matrices (dN, dS, omega)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(phangorn)
})

source(here("config.R"))
source(here("R", "matrix_utils.R"))
source(here("R", "tree_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] Assembling dN/dS matrices...")

  codeml_root <- file.path(PATHS$codeml, gene)
  acc_file    <- file.path(PATHS$metadata, paste0(gene, "_accessions_final.txt"))

  if (!dir.exists(codeml_root)) {
    warning("codeml directory not found, skipping: ", codeml_root)
    next
  }

  accessions <- readLines(acc_file)
  n <- length(accessions)
  message("  ", n, " accessions")

  # ── Collect pairwise matrices ─────────────────────────────────────────────
  pair_dirs <- list.dirs(codeml_root, full.names = TRUE, recursive = FALSE)
  pair_dirs <- pair_dirs[file.exists(file.path(pair_dirs, "2ML.dS"))]
  message("  ", length(pair_dirs), " completed codeml pairs found")

  dN_matrices <- lapply(pair_dirs, function(d) read_paml_matrix(file.path(d, "2ML.dN")))
  dS_matrices <- lapply(pair_dirs, function(d) read_paml_matrix(file.path(d, "2ML.dS")))

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

  # ── Build distance trees ──────────────────────────────────────────────────
  dN_tree <- nnls.tree(as.dist(dN_full), read.tree(file.path(PATHS$trees, paste0(gene, "_tree.nwk"))))
  dS_tree <- nnls.tree(as.dist(dS_full), read.tree(file.path(PATHS$trees, paste0(gene, "_tree.nwk"))))

  # ── Save outputs ──────────────────────────────────────────────────────────
  write.tree(dN_tree, file.path(PATHS$trees, paste0(gene, "_dN_tree.nwk")))
  write.tree(dS_tree, file.path(PATHS$trees, paste0(gene, "_dS_tree.nwk")))

  saveRDS(
    list(dN = dN_full, dS = dS_full, omega = dN_full / dS_full, accessions = common),
    file.path(PATHS$rates, paste0(gene, "_dNdS_matrix.rds"))
  )

  message("  dN tree:  ", file.path(PATHS$trees, paste0(gene, "_dN_tree.nwk")))
  message("  dS tree:  ", file.path(PATHS$trees, paste0(gene, "_dS_tree.nwk")))
  message("  Matrices: ", file.path(PATHS$rates,  paste0(gene, "_dNdS_matrix.rds")))
}

message("\nDone.")
