# scripts/01_alignment/02_macse.R
#
# PURPOSE: Run pairwise MACSE codon-aware alignments on accessions retained
#          after MAPLE tree pruning. Accessions are split into lettered batches
#          and all pairwise batch combinations are processed.
#
# INPUT:
#   data/metadata/{gene}_accessions_final.txt  — accessions after MAPLE pruning
#   data/raw/{gene}_sequences.fasta            — GISAID FASTA
#
# OUTPUT:
#   data/alignments/{gene}/macse/{pair}_aligned.rds  — list(dna_bin, pair_label)
#
# USAGE:
#   Rscript scripts/01_alignment/02_macse.R
#   Rscript scripts/01_alignment/02_macse.R --gene H1

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("config.R"))
source(here("R", "alignment_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] MACSE pairwise alignment...")

  acc_file <- file.path(PATHS$metadata, paste0(gene, "_accessions_final.txt"))
  fasta_in <- file.path(PATHS$raw,      paste0(gene, "_sequences.fasta"))
  out_dir  <- file.path(PATHS$alignments, gene, "macse")
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

  if (!file.exists(acc_file)) {
    warning("Accessions file not found — run 02_phylogenetics/01_maple_tree.R first: ", acc_file)
    next
  }

  accessions <- readLines(acc_file)
  message("  ", length(accessions), " accessions")

  # Split into batches
  batch_size <- PARAMS$macse_batch_size
  batch_labels <- LETTERS[seq_len(ceiling(length(accessions) / batch_size))]
  batches <- split(accessions, ceiling(seq_along(accessions) / batch_size))
  names(batches) <- batch_labels[seq_along(batches)]

  # All pairwise combinations (including self-pairs for within-batch)
  pairs <- combn(names(batches), 2, simplify = FALSE)
  # Add self-pairs
  self_pairs <- lapply(names(batches), function(x) c(x, x))
  pairs <- c(self_pairs, pairs)

  failed <- character(0)

  for (pair in pairs) {
    label <- paste0(pair[1], pair[2])
    out_file <- file.path(out_dir, paste0(label, "_aligned.rds"))
    if (file.exists(out_file)) {
      message("  [", label, "] Already done, skipping.")
      next
    }

    acc_subset <- unique(c(batches[[pair[1]]], batches[[pair[2]]]))

    result <- tryCatch(
      macse_align_pair(
        accessions = acc_subset,
        fasta_in   = fasta_in,
        out_dir    = out_dir,
        label      = label,
        macse_jar  = TOOLS$macse,
        mafft_bin  = TOOLS$mafft
      ),
      error = function(e) {
        message("  [", label, "] FAILED: ", conditionMessage(e))
        NULL
      }
    )

    if (!is.null(result)) {
      saveRDS(result, out_file)
      message("  [", label, "] Done.")
    } else {
      failed <- c(failed, label)
    }
  }

  if (length(failed) > 0) {
    fail_log <- file.path(out_dir, "failed_pairs.txt")
    writeLines(failed, fail_log)
    message("  ", length(failed), " failed pairs logged to: ", fail_log)
  }

  message("  Complete: ", length(pairs) - length(failed), "/", length(pairs), " pairs")
}

message("\nDone.")
