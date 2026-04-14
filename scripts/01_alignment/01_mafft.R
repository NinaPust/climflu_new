# scripts/01_alignment/01_mafft.R
#
# PURPOSE: Run MAFFT multiple sequence alignment on each gene segment,
#          then trim codon columns with excessive gaps.
#
# INPUT:
#   data/raw/{gene}_sequences.fasta  — GISAID FASTA export
#
# OUTPUT:
#   data/alignments/{gene}_aligned.fasta  — trimmed alignment

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("config.R"))
source(here("R", "alignment_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] MAFFT alignment...")

  fasta_in  <- file.path(PATHS$raw,        paste0(gene, "_sequences.fasta"))
  fasta_out <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))

  if (!file.exists(fasta_in)) {
    warning("FASTA not found, skipping: ", fasta_in)
    next
  }

  aligned <- run_mafft(fasta_in, fasta_out, threads = PARAMS$mafft_threads,
                       mafft_bin = TOOLS$mafft)

  trimmed <- trim_codons(aligned, gap_threshold = PARAMS$trim_gap_threshold)

  write.FASTA(trimmed, fasta_out)
  message("  Written: ", fasta_out, " (", length(trimmed), " sequences, ",
          ncol(as.matrix(trimmed)), " sites)")
}

message("\nDone.")
