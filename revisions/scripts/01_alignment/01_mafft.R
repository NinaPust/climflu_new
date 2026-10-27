# revisions/scripts/01_alignment/01_mafft.R
#
# PURPOSE (REVISION MODE — differs from the original pipeline script):
#   Reproduce the PUBLISHED alignment and quantify how many codon columns the
#   >50% gap trim discards (Reviewer 1, major comment 2).
#
# WHY THIS DIFFERS FROM THE ORIGINAL:
#   The original script filters the raw FASTA to accessions present in
#   data/climate/{gene}_climate.csv. That climate CSV was regenerated on
#   Apr 27, AFTER the published alignments were built on Apr 16, so the filter
#   no longer selects the same sequence set. Re-running the original logic today
#   yields a DIFFERENT alignment (for H1: 10,612 sequences instead of 10,918,
#   with different columns trimmed and every row shifted), which would not
#   correspond to the existing MAPLE trees, TreeTime chronograms, or rates.
#
#   The published accession set is still recoverable: it is exactly the set of
#   sequence names in the published alignments themselves. This
#   script therefore takes its sequence set from there rather than from the
#   climate CSV, so the re-run is anchored to what was actually published.
#
# SAFETY:
#   - Reads the published alignments read-only; never writes to data/alignments/.
#   - Writes to revisions/data/alignments_repro/ — never to data/alignments/,
#     so a re-run cannot clobber either the published or the current alignments.
#   - Diffs the reproduced alignment against the published one and reports
#     whether they match, so you know whether the trim counts apply exactly.
#
# NOTE ON DEDUPLICATION:
#   The published set has already passed both the pre- and post-alignment dedup
#   steps, so this script does NOT dedup again — doing so could drop sequences
#   and break the 1:1 correspondence with the published set. Any duplicates that
#   do appear are reported as a diagnostic instead.
#
# INPUT:
#   data/raw/{gene}_sequences.fasta
#   data/alignments/{gene}_aligned.fasta          (published, read-only)
#
# OUTPUT:
#   revisions/data/alignments_repro/{gene}_aligned_repro.fasta
#   revisions/data/alignments_repro/trim_stats.csv
#
# FLAGS:
#   --gene X   run a single gene (default: all six)
#   --force    overwrite an existing reproduced alignment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "alignment_utils.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else ALL_GENES
force <- "--force" %in% args

# The published alignments live in PATHS$alignments. (They were temporarily
# relocated to alignments/old/ during the MAFFT re-run; that is no longer the
# case.) Override with --published-dir if they are ever moved again.
published_dir <- if ("--published-dir" %in% args)
  args[which(args == "--published-dir") + 1] else PATHS$alignments
out_dir       <- here("revisions", "data", "alignments_repro")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stats_rows <- list()

for (gene in genes) {
  message("\n[", gene, "] ═══════════════════════════════════════")

  fasta_in   <- file.path(PATHS$raw, paste0(gene, "_sequences.fasta"))
  published  <- file.path(published_dir, paste0(gene, "_aligned.fasta"))
  fasta_out  <- file.path(out_dir, paste0(gene, "_aligned_repro.fasta"))

  if (file.exists(fasta_out) && !force) {
    message("  Reproduced alignment exists, skipping (use --force): ", fasta_out)
    next
  }
  if (!file.exists(fasta_in))  { warning("Raw FASTA not found: ", fasta_in);   next }
  if (!file.exists(published)) { warning("Published alignment not found: ",
                                         published); next }

  # ── Step 1: recover the published accession set ──────────────────────────
  pub_aln <- read.FASTA(published)
  pub_acc <- names(pub_aln)
  pub_mat <- as.matrix(read.dna(published, format = "fasta",
                                as.character = TRUE, as.matrix = TRUE))
  message("  Published alignment: ", length(pub_acc), " sequences, ",
          ncol(pub_mat), " sites")

  # ── Step 2: subset the raw FASTA to exactly those accessions ─────────────
  all_seqs <- read.FASTA(fasta_in)
  message("  Raw FASTA: ", length(all_seqs), " sequences")
  missing <- setdiff(pub_acc, names(all_seqs))
  if (length(missing) > 0)
    warning("  ", length(missing), " published accessions absent from raw FASTA")
  seqs <- all_seqs[intersect(pub_acc, names(all_seqs))]
  message("  Subset to published set: ", length(seqs), " sequences")
  if (length(seqs) < 2) { warning("Too few sequences for ", gene); next }

  # ── Step 3: MAFFT (no dedup — see header) ────────────────────────────────
  filtered_fasta <- file.path(out_dir, paste0(gene, "_reproinput.fasta"))
  write.FASTA(seqs, filtered_fasta)

  message("  Running MAFFT...")
  t0 <- proc.time()
  raw_out <- file.path(out_dir, paste0(gene, "_untrimmed.fasta"))
  aligned <- run_mafft(filtered_fasta, raw_out,
                       threads   = PARAMS$mafft_threads,
                       mafft_bin = TOOLS$mafft)
  message("  MAFFT done in ", round((proc.time() - t0)[["elapsed"]] / 60, 1),
          " min")

  # ── Step 4: trim gappy codon columns, capturing the counts ───────────────
  trimmed <- trim_codons(aligned, gap_threshold = PARAMS$trim_gap_threshold)
  ts <- attr(trimmed, "trim_stats")
  message("  TRIM: ", ts$codons_removed, " of ", ts$codons_before,
          " codons removed (", ts$pct_codons_removed, "%), ",
          ts$codons_after, " retained")

  trimmed_mat <- as.matrix(trimmed)

  n_dup <- sum(duplicated(apply(trimmed_mat, 1, paste, collapse = "")))
  if (n_dup > 0)
    message("  NOTE: ", n_dup, " duplicate rows present after trimming ",
            "(not removed — see header)")

  write.FASTA(as.list(trimmed_mat), fasta_out)
  message("  Written: ", fasta_out, " (", nrow(trimmed_mat), " sequences, ",
          ncol(trimmed_mat), " sites)")

  # ── Step 5: verify against the published alignment ───────────────────────
  same_dim  <- nrow(trimmed_mat) == nrow(pub_mat) &&
               ncol(trimmed_mat) == ncol(pub_mat)
  shared    <- intersect(rownames(trimmed_mat), rownames(pub_mat))
  identical_rows <- if (ncol(trimmed_mat) == ncol(pub_mat)) {
    sum(vapply(shared, function(s)
      identical(toupper(unname(trimmed_mat[s, ])),
                toupper(unname(pub_mat[s, ]))), logical(1)))
  } else NA_integer_

  message("  VERIFY vs published: dims ", if (same_dim) "MATCH" else "DIFFER",
          " (", nrow(trimmed_mat), "x", ncol(trimmed_mat), " vs ",
          nrow(pub_mat), "x", ncol(pub_mat), "); ",
          "identical rows: ", identical_rows, "/", length(shared))
  if (isTRUE(same_dim) && identical(identical_rows, length(shared))) {
    message("  => REPRODUCED EXACTLY. Trim counts apply to the published ",
            "alignment.")
  } else {
    message("  => NOT an exact reproduction. Report trim counts as ",
            "approximate for the published alignment.")
  }

  stats_rows[[gene]] <- data.frame(
    gene                = gene,
    published_seqs      = nrow(pub_mat),
    published_sites     = ncol(pub_mat),
    repro_seqs          = nrow(trimmed_mat),
    repro_sites         = ncol(trimmed_mat),
    codons_before_trim  = ts$codons_before,
    codons_removed      = ts$codons_removed,
    codons_after_trim   = ts$codons_after,
    pct_codons_removed  = ts$pct_codons_removed,
    nt_before_trim      = ts$nt_before,
    nt_after_trim       = ts$nt_after,
    gap_threshold       = ts$gap_threshold,
    shared_seqs         = length(shared),
    identical_rows      = identical_rows,
    exact_reproduction  = isTRUE(same_dim) &&
                          identical(identical_rows, length(shared)),
    stringsAsFactors    = FALSE
  )

  unlink(filtered_fasta)
}

if (length(stats_rows) > 0) {
  stats <- do.call(rbind, stats_rows)
  stats_path <- file.path(out_dir, "trim_stats.csv")
  if (file.exists(stats_path)) {
    prev <- read_csv(stats_path, show_col_types = FALSE)
    stats <- rbind(prev[!prev$gene %in% stats$gene, , drop = FALSE], stats)
  }
  write_csv(stats, stats_path)
  cat("\n================ TRIM SUMMARY ================\n")
  print(as.data.frame(stats[, c("gene", "codons_before_trim", "codons_removed",
                                "pct_codons_removed", "exact_reproduction")]),
        row.names = FALSE)
  message("\nWritten: ", stats_path)
}

message("\nDone.")
