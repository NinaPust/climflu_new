# R/alignment_utils.R
# Utilities for sequence alignment (MAFFT) and codon trimming (MACSE).

#' Run MAFFT on a FASTA file and return the aligned DNAbin object.
#'
#' @param fasta_in   Path to input FASTA.
#' @param fasta_out  Path to write aligned FASTA.
#' @param threads    Number of CPU threads for MAFFT.
#' @param mafft_bin  Path to mafft executable.
#' @return DNAbin object of the alignment.
run_mafft <- function(fasta_in, fasta_out, threads = 4, mafft_bin = "mafft") {
  if (!file.exists(fasta_in)) stop("Input FASTA not found: ", fasta_in)

  cmd_args <- c("--auto", "--thread", threads, "--out", fasta_out, fasta_in)
  ret <- system2(mafft_bin, args = cmd_args, stdout = TRUE, stderr = TRUE)

  if (!file.exists(fasta_out)) stop("MAFFT did not produce output: ", fasta_out)

  ape::read.FASTA(fasta_out)
}


#' Trim alignment columns where gap frequency exceeds a threshold.
#'
#' Operates at the codon level: only trims complete codons (3 columns),
#' removing a codon if any of its three positions exceeds the gap threshold.
#'
#' @param alignment     DNAbin alignment matrix.
#' @param gap_threshold Proportion of gaps above which a codon is removed (default 0.5).
#' @return Trimmed DNAbin alignment.
trim_codons <- function(alignment, gap_threshold = 0.5) {
  mat <- as.matrix(alignment)
  n_seq  <- nrow(mat)
  n_cols <- ncol(mat)

  # Convert DNAbin matrix to character matrix for gap detection
  char_mat <- as.character(mat)
  gap_freq <- colMeans(char_mat == "-", na.rm = TRUE)

  # Mark codons to keep: all 3 positions must be below threshold
  n_codons   <- floor(n_cols / 3)
  keep_cols  <- logical(n_cols)

  for (k in seq_len(n_codons)) {
    cols <- ((k - 1) * 3 + 1):(k * 3)
    if (all(cols <= n_cols) && all(gap_freq[cols] <= gap_threshold)) {
      keep_cols[cols] <- TRUE
    }
  }

  n_removed <- sum(!keep_cols)
  if (n_removed > 0) {
    message("  Trimmed ", n_removed / 3, " codon columns (", n_removed, " nt) ",
            "with gap frequency > ", gap_threshold)
  }

  ape::as.DNAbin(char_mat[, keep_cols, drop = FALSE])
}


#' Run MACSE codon-aware alignment on a subset of sequences.
#'
#' Writes a temporary FASTA, aligns with MAFFT, then refines with MACSE.
#'
#' @param accessions  Character vector of accession IDs.
#' @param fasta_in    Path to full FASTA file to extract sequences from.
#' @param out_dir     Directory to write temporary files.
#' @param label       Label for this pair (used in filenames).
#' @param macse_jar   Path to MACSE JAR file.
#' @param mafft_bin   Path to mafft executable.
#' @return DNAbin object of MACSE-aligned sequences.
macse_align_pair <- function(accessions, fasta_in, out_dir, label,
                              macse_jar, mafft_bin = "mafft") {
  all_seqs <- ape::read.FASTA(fasta_in)
  present  <- intersect(accessions, names(all_seqs))
  if (length(present) < 2) stop("Fewer than 2 sequences found for pair ", label)

  subset_seqs <- all_seqs[present]

  # Write subset FASTA
  tmp_fasta  <- file.path(out_dir, paste0(label, "_input.fasta"))
  mafft_out  <- file.path(out_dir, paste0(label, "_mafft.fasta"))
  macse_out  <- file.path(out_dir, paste0(label, "_macse"))

  ape::write.FASTA(subset_seqs, tmp_fasta)

  # MAFFT pre-alignment
  run_mafft(tmp_fasta, mafft_out, threads = 1, mafft_bin = mafft_bin)

  # MACSE refinement
  macse_args <- c(
    "-jar",  macse_jar,
    "-prog", "refineAlignment",
    "-align", mafft_out,
    "-out_NT", paste0(macse_out, "_NT.fasta"),
    "-out_AA", paste0(macse_out, "_AA.fasta")
  )
  ret <- system2("java", args = macse_args, stdout = TRUE, stderr = TRUE)

  nt_file <- paste0(macse_out, "_NT.fasta")
  if (!file.exists(nt_file)) stop("MACSE did not produce NT output for: ", label)

  ape::read.FASTA(nt_file)
}
