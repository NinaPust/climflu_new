# scripts/04_epitope_analysis/02_split_alignments.R
#
# PURPOSE: Split the MAFFT HA alignment into epitope, non-epitope, ridge,
#          and non-ridge sub-alignments using column indices from
#          01_map_epitope_sites.R. Only sequences present in the post-MAPLE
#          pruned accession list are retained.
#
# INPUT:
#   data/alignments/{gene}_aligned.fasta
#   data/trees/pruned/{gene}_accessions.txt
#   data/epitope/{gene}_epitope_cols.txt
#   data/epitope/{gene}_nonepitope_cols.txt
#   data/epitope/{gene}_ridge_cols.txt
#   data/epitope/{gene}_nonridge_cols.txt
#
# OUTPUT:
#   data/epitope/{gene}_epitope_aligned.fasta
#   data/epitope/{gene}_nonepitope_aligned.fasta
#   data/epitope/{gene}_ridge_aligned.fasta
#   data/epitope/{gene}_nonridge_aligned.fasta
#
# FLAGS:
#   --gene X   run a single gene (default: H1 H3)

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("config.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  c("H1", "H3")
}

epi_dir <- here("data", "epitope")

for (gene in genes) {
  message("\n[", gene, "] Splitting alignment...")

  aln_path     <- file.path(PATHS$alignments,  paste0(gene, "_aligned.fasta"))
  acc_path     <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
  epi_col_path   <- file.path(epi_dir, paste0(gene, "_epitope_cols.txt"))
  non_col_path   <- file.path(epi_dir, paste0(gene, "_nonepitope_cols.txt"))
  ridge_col_path <- file.path(epi_dir, paste0(gene, "_ridge_cols.txt"))
  nridge_col_path <- file.path(epi_dir, paste0(gene, "_nonridge_cols.txt"))

  for (p in c(aln_path, acc_path, epi_col_path, non_col_path,
              ridge_col_path, nridge_col_path)) {
    if (!file.exists(p)) stop("Required file not found: ", p)
  }

  # Load alignment as character matrix for column subsetting
  message("  Reading alignment...")
  aln <- read.dna(aln_path, format = "fasta",
                  as.character = TRUE, as.matrix = TRUE)
  message("  ", nrow(aln), " sequences, ", ncol(aln), " columns")

  # Subset to pruned accessions
  accessions <- readLines(acc_path)
  accessions <- accessions[nchar(trimws(accessions)) > 0]
  common     <- intersect(accessions, rownames(aln))
  message("  ", length(common), " sequences after pruning accession filter (",
          length(accessions) - length(common), " not in alignment)")
  aln <- aln[common, , drop = FALSE]

  # Read column indices (1-based nt columns)
  epi_cols   <- as.integer(readLines(epi_col_path))
  non_cols   <- as.integer(readLines(non_col_path))
  ridge_cols <- as.integer(readLines(ridge_col_path))
  nridge_cols <- as.integer(readLines(nridge_col_path))
  message("  Epitope columns: ", length(epi_cols),
          "  Non-epitope columns: ", length(non_cols),
          "  Ridge columns: ", length(ridge_cols))

  # Subset and remove all-gap sequences
  split_and_write <- function(cols, label, out_path) {
    sub_aln <- aln[, cols, drop = FALSE]
    gap_frac <- rowMeans(sub_aln == "-" | sub_aln == "n" | sub_aln == "?")
    all_gap  <- which(gap_frac == 1)
    if (length(all_gap) > 0) {
      message("  [", label, "] Dropping ", length(all_gap),
              " all-gap sequences")
      sub_aln <- sub_aln[-all_gap, , drop = FALSE]
    }
    message("  [", label, "] Writing ", nrow(sub_aln),
            " sequences, ", ncol(sub_aln), " columns → ", basename(out_path))
    write.dna(as.DNAbin(sub_aln), out_path, format = "fasta",
              colsep = "", nbcol = -1)
  }

  split_and_write(epi_cols, "epitope",
                  file.path(epi_dir, paste0(gene, "_epitope_aligned.fasta")))
  split_and_write(non_cols, "nonepitope",
                  file.path(epi_dir, paste0(gene, "_nonepitope_aligned.fasta")))
  split_and_write(ridge_cols, "ridge",
                  file.path(epi_dir, paste0(gene, "_ridge_aligned.fasta")))
  split_and_write(nridge_cols, "nonridge",
                  file.path(epi_dir, paste0(gene, "_nonridge_aligned.fasta")))
}

message("\nDone.")
