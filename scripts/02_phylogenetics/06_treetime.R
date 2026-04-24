# scripts/02_phylogenetics/06_treetime.R
#
# PURPOSE: Date the canonical (r2) pruned MAPLE tree using TreeTime.
#          This is the timetree used for all downstream analysis.
#
# INPUT:
#   data/trees/pruned/{gene}_pruned.nwk          — canonical r2 pruned tree
#   data/alignments/{gene}_aligned.fasta
#   data/metadata/{gene}_metadata.csv
#
# OUTPUT:
#   data/trees/timed/{gene}_timetree.nwk         — time-calibrated tree
#
# REQUIRES: TOOLS$treetime in config.R
#           (treetime_env: /home/nina/miniforge3/envs/treetime_env/bin/treetime)
#
# FLAGS:
#   --gene X   run a single gene segment
#   --test     use subsampled data from data/test/ (see 00_subsample_test.R)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
  library(dplyr)
})

source(here("config.R"))

args     <- commandArgs(trailingOnly = TRUE)
ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
genes    <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}
use_test <- "--test" %in% args

if (use_test) message("*** Running in TEST mode (subsampled data) ***")

paths_pruned <- if (use_test) PATHS$test_trees_pruned else PATHS$trees_pruned
paths_meta   <- if (use_test) PATHS$test_metadata     else PATHS$metadata
paths_timed  <- if (use_test) PATHS$test_trees_timed  else PATHS$trees_timed

for (gene in genes) {
  message("\n[", gene, "] TreeTime molecular clock dating...")

  tree_in  <- file.path(paths_pruned,    paste0(gene, "_pruned.nwk"))
  aln_in   <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  meta_in  <- file.path(paths_meta,      paste0(gene, "_metadata.csv"))
  tree_out <- file.path(paths_timed,     paste0(gene, "_timetree.nwk"))
  tt_dir   <- file.path(paths_timed,     paste0(gene, "_treetime_out"))

  dir.create(paths_timed, recursive = TRUE, showWarnings = FALSE)

  if (!file.exists(tree_in)) {
    warning("Tree not found: ", tree_in)
    next
  }
  if (!file.exists(aln_in)) {
    warning("Alignment not found: ", aln_in)
    next
  }
  if (!file.exists(meta_in)) {
    warning("Metadata not found: ", meta_in)
    next
  }

  dir.create(tt_dir, recursive = TRUE, showWarnings = FALSE)

  # Write dates file (CSV: name,date) — treetime accepts ISO dates
  meta  <- read_csv(meta_in, show_col_types = FALSE)
  dates <- meta |>
    select(name = accession, date = collection_date) |>
    filter(!is.na(date)) |>
    mutate(date = format(as.Date(date), "%Y-%m-%d"))

  date_file <- file.path(tt_dir, "dates.csv")
  write_csv(dates, date_file)

  cmd <- paste(
    shQuote(TOOLS$treetime),
    "--tree",   shQuote(tree_in),
    "--aln",    shQuote(aln_in),
    "--dates",  shQuote(date_file),
    "--relax",  "1.0 0.0",
    "--outdir", shQuote(tt_dir)
  )
  message("  CMD: ", cmd)
  exit_code <- system(cmd)

  if (exit_code != 0) {
    warning("TreeTime returned non-zero exit code for ", gene,
            " (exit ", exit_code, ")")
    next
  }

  # TreeTime writes timetree.nexus (or .nwk) depending on version
  tt_files <- list.files(tt_dir,
                         pattern    = "timetree\\.(nwk|nexus|tre)$",
                         full.names = TRUE)
  if (length(tt_files) == 0) {
    warning("TreeTime produced no timetree output for ", gene,
            ". Check: ", tt_dir)
    next
  }

  tt_file  <- tt_files[1]
  timetree <- if (grepl("\\.nexus$", tt_file)) {
    read.nexus(tt_file)
  } else {
    read.tree(tt_file)
  }

  # Remove clock outliers flagged by TreeTime
  outlier_file <- file.path(tt_dir, "outliers.tsv")
  if (file.exists(outlier_file)) {
    outliers <- read.delim(outlier_file, stringsAsFactors = FALSE)
    bad_tips <- outliers[[1]][!is.na(outliers[[1]]) & outliers[[1]] != ""]
    bad_tips <- intersect(bad_tips, timetree$tip.label)
    if (length(bad_tips) > 0) {
      message("  Removing ", length(bad_tips), " clock outlier(s): ",
              paste(bad_tips, collapse = ", "))
      timetree <- drop.tip(timetree, bad_tips)
    } else {
      message("  No clock outliers to remove.")
    }
  }

  write.tree(timetree, tree_out)
  message("  Written: ", tree_out, " (", Ntip(timetree), " tips)")
}

message("\nDone.")
