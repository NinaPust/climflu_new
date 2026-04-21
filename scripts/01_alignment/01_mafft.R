# scripts/01_alignment/01_mafft.R
#
# PURPOSE: Filter sequences to those with climate data, remove duplicate
#          sequences, align with MAFFT, trim gappy codons, remove
#          post-alignment duplicates, and sync the climate CSV.
#
# PIPELINE:
#   1. Load climate CSV → get accessions with climate data
#   2. Subset raw FASTA to those accessions
#   3. Remove duplicate sequences (keep first occurrence)
#   4. MAFFT alignment
#   5. Trim codon columns with excessive gaps
#   6. Remove duplicate sequences again (alignment may reveal identity)
#
# NOTE: Climate CSVs are NOT modified here. Genes within the same subtype
#       (H1/N1/PA_H1N1 and H3/N2/PA_H3N2) share metadata, so dropping
#       accessions per-gene would desync them.
#
# INPUT:
#   data/raw/{gene}_sequences.fasta
#   data/climate/{gene}_climate.csv
#
# OUTPUT:
#   data/alignments/{gene}_aligned.fasta  — final trimmed alignment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
})

source(here("config.R"))
source(here("R", "alignment_utils.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else ALL_GENES

#' Remove duplicate sequences from a DNAbin list.
#' Keeps the first occurrence of each unique sequence.
#' @param seqs DNAbin list
#' @return DNAbin list with duplicates removed
remove_duplicate_seqs <- function(seqs) {
  # Convert each sequence to a character string for comparison
  seq_strings <- vapply(seqs, function(s) {
    paste(as.character(s), collapse = "")
  }, character(1))

  dup <- duplicated(seq_strings)
  n_dup <- sum(dup)
  if (n_dup > 0) {
    message("  Removed ", n_dup, " duplicate sequences (",
            length(seqs) - n_dup, " unique)")
  }
  seqs[!dup]
}

for (gene in genes) {
  message("\n[", gene, "] ═══════════════════════════════════════")

  fasta_in    <- file.path(PATHS$raw, paste0(gene, "_sequences.fasta"))
  climate_csv <- file.path(PATHS$climate, paste0(gene, "_climate.csv"))
  fasta_out   <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))

  if (file.exists(fasta_out)) {
    message("  Alignment already exists, skipping: ", fasta_out)
    next
  }
  if (!file.exists(fasta_in)) {
    warning("FASTA not found, skipping: ", fasta_in)
    next
  }
  if (!file.exists(climate_csv)) {
    warning("Climate CSV not found, skipping: ", climate_csv)
    next
  }

  # ── Step 1: get accessions with climate data ──────────────────────────────
  climate <- read_csv(climate_csv, show_col_types = FALSE)
  climate_acc <- climate$accession
  message("  Climate CSV: ", length(climate_acc), " accessions")

  # ── Step 2: subset FASTA to climate accessions ────────────────────────────
  all_seqs <- read.FASTA(fasta_in)
  message("  Raw FASTA: ", length(all_seqs), " sequences")

  keep <- intersect(names(all_seqs), climate_acc)
  seqs <- all_seqs[keep]
  message("  After climate filter: ", length(seqs), " sequences")

  if (length(seqs) < 2) {
    warning("Fewer than 2 sequences after filtering for ", gene)
    next
  }

  # ── Step 3: remove duplicate sequences (pre-alignment) ───────────────────
  message("  Pre-alignment dedup:")
  seqs <- remove_duplicate_seqs(seqs)

  # Write filtered FASTA for MAFFT input
  filtered_fasta <- file.path(PATHS$alignments,
                              paste0(gene, "_filtered.fasta"))
  write.FASTA(seqs, filtered_fasta)

  # ── Step 4: MAFFT alignment ──────────────────────────────────────────────
  message("  Running MAFFT...")
  aligned <- run_mafft(filtered_fasta, fasta_out,
                       threads = PARAMS$mafft_threads,
                       mafft_bin = TOOLS$mafft)

  # ── Step 5: trim gappy codon columns ─────────────────────────────────────
  trimmed <- trim_codons(aligned,
                         gap_threshold = PARAMS$trim_gap_threshold)

  # ── Step 6: remove duplicate sequences (post-alignment) ──────────────────
  message("  Post-alignment dedup:")
  trimmed_list <- remove_duplicate_seqs(as.list(trimmed))
  # Rebuild aligned matrix from kept sequences
  trimmed_mat <- as.matrix(trimmed)[names(trimmed_list), , drop = FALSE]

  write.FASTA(as.list(trimmed_mat), fasta_out)
  message("  Alignment written: ", fasta_out, " (",
          nrow(trimmed_mat), " sequences, ",
          ncol(trimmed_mat), " sites)")

  # Clean up intermediate file
  unlink(filtered_fasta)
}

message("\nDone.")
