# scripts/04_epitope_analysis/04_assemble_matrices.R
#
# PURPOSE: Stitch per-pair yn00 RDS files into full symmetric dN/dS
#          matrices for epitope and non-epitope sub-alignments.
#          Pair directories are deleted after stitching to free disk space.
#          Mirrors 02_phylogenetics/05_assemble_matrices.R.
#
# INPUT:   data/epitope/yn00/{gene}_{type}/*/dN_remapped.rds
# OUTPUT:  data/epitope/{gene}_{type}_dNdS_matrix.rds
#
# FLAGS:
#   --gene X   run a single gene (default: H1 H3)
#   --type X   epitope|nonepitope (default: both)

suppressPackageStartupMessages({
  library(here)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "matrix_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  c("H1", "H3")
}
types <- if ("--type" %in% args) {
  args[which(args == "--type") + 1]
} else {
  c("epitope", "nonepitope", "ridge", "nonridge")
}

epi_dir <- here("revisions", "data", "epitope")

for (gene in genes) {
  for (type in types) {
    tag <- paste0(gene, "_", type)
    message("\n[", tag, "] Assembling dN/dS matrices...")

    yn00_root <- file.path(epi_dir, "yn00", tag)
    if (!dir.exists(yn00_root)) {
      warning("yn00 directory not found, skipping: ", yn00_root)
      next
    }

    # Accession list: all sequences present in the sub-alignment FASTA
    fasta_path <- file.path(epi_dir, paste0(tag, "_aligned.fasta"))
    if (!file.exists(fasta_path)) {
      warning("Sub-alignment FASTA not found: ", fasta_path)
      next
    }
    aln        <- ape::read.FASTA(fasta_path)
    accessions <- names(aln)
    rm(aln)
    message("  ", length(accessions), " accessions from sub-alignment")

    # Collect completed pair directories (have at least dS_remapped.rds)
    pair_dirs <- list.dirs(yn00_root, full.names = TRUE, recursive = FALSE)
    pair_dirs <- pair_dirs[
      file.exists(file.path(pair_dirs, "dS_remapped.rds"))
    ]
    message("  ", length(pair_dirs), " completed pairs found")

    if (length(pair_dirs) == 0) {
      warning("No completed pairs found for ", tag, " — skipping.")
      next
    }

    load_rds <- function(d, type_str) {
      path <- file.path(d, paste0(type_str, "_remapped.rds"))
      tryCatch(readRDS(path),
               error = function(e) {
                 warning("Could not read ", path, ": ", e$message)
                 NULL
               })
    }

    dn_matrices <- lapply(pair_dirs, load_rds, type_str = "dN")
    ds_matrices <- lapply(pair_dirs, load_rds, type_str = "dS")

    dn_full <- matrix_stitch(dn_matrices, expected_labels = accessions)
    ds_full <- matrix_stitch(ds_matrices, expected_labels = accessions)

    common  <- intersect(rownames(dn_full), accessions)
    removed <- setdiff(accessions, common)
    if (length(removed) > 0)
      message("  ", length(removed),
              " accessions absent from all pairs (all NA rows)")

    dn_full <- dn_full[common, common]
    ds_full <- ds_full[common, common]

    validate_matrix(dn_full, expected_n = length(common))
    validate_matrix(ds_full, expected_n = length(common))

    out_path <- file.path(epi_dir, paste0(tag, "_dNdS_matrix.rds"))
    saveRDS(
      list(dN         = dn_full,
           dS         = ds_full,
           omega      = dn_full / pmax(ds_full, 1e-10),
           accessions = common),
      out_path
    )
    message("  Saved: ", out_path)

    # Delete pair directories now that matrices are stitched
    message("  Deleting ", length(pair_dirs), " pair directories...")
    unlink(pair_dirs, recursive = TRUE)
    message("  Done.")
  }
}

message("\nDone.")
