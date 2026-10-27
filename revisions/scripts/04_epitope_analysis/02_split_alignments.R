# scripts/04_epitope_analysis/02_split_alignments.R
#
# PURPOSE: Split the MAFFT HA alignment into epitope, non-epitope, ridge,
#          non-ridge, and HA2 sub-alignments using column indices from
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
#   data/epitope/{gene}_site_map.csv              (for HA2 column derivation)
#
# OUTPUT:
#   data/epitope/{gene}_epitope_aligned.fasta
#   data/epitope/{gene}_nonepitope_aligned.fasta
#   data/epitope/{gene}_ridge_aligned.fasta
#   data/epitope/{gene}_nonridge_aligned.fasta
#   data/epitope/{gene}_HA2_aligned.fasta
#   data/epitope/{gene}_HA2_cols.txt
#
# FLAGS:
#   --gene X    run a single gene (default: H1 H3)
#   --types X   comma-separated subset of types to run (default: all)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
  library(dplyr)
})

source(here("revisions", "config.R"))

args       <- commandArgs(trailingOnly = TRUE)
genes      <- if ("--gene"  %in% args) args[which(args == "--gene")  + 1] else c("H1", "H3")
types_arg  <- if ("--types" %in% args) strsplit(args[which(args == "--types") + 1], ",")[[1]] else NULL
all_types  <- c("epitope", "nonepitope", "ridge", "nonridge", "HA2", "full")
run_types  <- if (!is.null(types_arg)) intersect(types_arg, all_types) else all_types

epi_dir <- here("revisions", "data", "epitope")

for (gene in genes) {
  message("\n[", gene, "] Splitting alignment...")

  aln_path     <- file.path(PATHS$alignments,  paste0(gene, "_aligned.fasta"))
  acc_path     <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
  epi_col_path    <- file.path(epi_dir, paste0(gene, "_epitope_cols.txt"))
  non_col_path    <- file.path(epi_dir, paste0(gene, "_nonepitope_cols.txt"))
  ridge_col_path  <- file.path(epi_dir, paste0(gene, "_ridge_cols.txt"))
  nridge_col_path <- file.path(epi_dir, paste0(gene, "_nonridge_cols.txt"))
  site_map_path   <- file.path(epi_dir, paste0(gene, "_site_map.csv"))

  required <- c(aln_path, acc_path)
  if (any(c("epitope", "nonepitope") %in% run_types))
    required <- c(required, epi_col_path, non_col_path)
  if (any(c("ridge", "nonridge") %in% run_types))
    required <- c(required, ridge_col_path, nridge_col_path)
  if (any(c("HA2", "full") %in% run_types))
    required <- c(required, site_map_path)
  for (p in required) {
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

  if ("epitope" %in% run_types || "nonepitope" %in% run_types) {
    epi_cols <- as.integer(readLines(epi_col_path))
    non_cols <- as.integer(readLines(non_col_path))
    message("  Epitope columns: ", length(epi_cols),
            "  Non-epitope columns: ", length(non_cols))
    if ("epitope" %in% run_types)
      split_and_write(epi_cols, "epitope",
                      file.path(epi_dir, paste0(gene, "_epitope_aligned.fasta")))
    if ("nonepitope" %in% run_types)
      split_and_write(non_cols, "nonepitope",
                      file.path(epi_dir, paste0(gene, "_nonepitope_aligned.fasta")))
  }

  if ("ridge" %in% run_types || "nonridge" %in% run_types) {
    ridge_cols  <- as.integer(readLines(ridge_col_path))
    nridge_cols <- as.integer(readLines(nridge_col_path))
    message("  Ridge columns: ", length(ridge_cols))
    if ("ridge" %in% run_types)
      split_and_write(ridge_cols, "ridge",
                      file.path(epi_dir, paste0(gene, "_ridge_aligned.fasta")))
    if ("nonridge" %in% run_types)
      split_and_write(nridge_cols, "nonridge",
                      file.path(epi_dir, paste0(gene, "_nonridge_aligned.fasta")))
  }

  if (any(c("HA2", "full") %in% run_types))
    sm <- read_csv(site_map_path, show_col_types = FALSE)

  if ("HA2" %in% run_types) {
    ha2_cols <- sm |> filter(region == "HA2") |> pull(alignment_col) |> sort()
    message("  HA2 columns (nt): ", length(ha2_cols))
    cols_out <- file.path(epi_dir, paste0(gene, "_HA2_cols.txt"))
    writeLines(as.character(ha2_cols), cols_out)
    split_and_write(ha2_cols, "HA2",
                    file.path(epi_dir, paste0(gene, "_HA2_aligned.fasta")))
  }

  # ── full: the complete CDS, in frame ───────────────────────────────────────
  # Previously 09_substitutions.R read the RAW whole-gene alignment for this
  # region and split it into triplets from column 1. That alignment begins with
  # a non-coding 5' leader (and ends with a trailer), so every `full` codon was
  # out of frame — the same defect that affected the other regions, but via a
  # different route, so the 01_map_epitope_sites.R fix did not reach it.
  #
  # Taking the columns from the site map instead gives exactly the reference's
  # CDS, in codon order, with the non-coding flanks excluded. Note this changes
  # what `full` means relative to the published analysis, which treated the
  # entire alignment (leader and trailer included) as coding sequence.
  #
  # Definition here is the COMPLETE CDS: signal peptide + HA1 + HA2. To restrict
  # to the mature protein instead, filter out region == "signal" below.
  if ("full" %in% run_types) {
    full_cols <- sm |> pull(alignment_col) |> unique() |> sort()
    n_codons  <- length(full_cols) / 3
    message("  Full CDS columns (nt): ", length(full_cols),
            "  (", n_codons, " codons)")
    if (length(full_cols) %% 3 != 0)
      warning("full CDS column count (", length(full_cols),
              ") is not a multiple of 3 — codon grouping will be wrong")
    cols_out <- file.path(epi_dir, paste0(gene, "_full_cols.txt"))
    writeLines(as.character(full_cols), cols_out)
    split_and_write(full_cols, "full",
                    file.path(epi_dir, paste0(gene, "_full_aligned.fasta")))
  }
}

message("\nDone.")
