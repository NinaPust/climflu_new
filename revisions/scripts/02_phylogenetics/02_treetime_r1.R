# scripts/02_phylogenetics/02_treetime_r1.R
#
# PURPOSE: Run TreeTime on the r1 pruned tree to identify clock outliers.
#          The output (outliers.tsv) is consumed by 03_maple_r2.R to exclude
#          temporally misplaced sequences before rebuilding the tree.
#          The r1 timetree is NOT used for downstream analysis.
#
# INPUT:
#   data/trees/pruned/{gene}_pruned_r1.nwk
#   data/alignments/{gene}_aligned.fasta
#   data/metadata/{gene}_metadata.csv
#
# OUTPUT:
#   data/trees/timed/{gene}_timetree_r1.nwk       — r1 timetree (QC only)
#   data/trees/timed/{gene}_treetime_r1_out/       — full TreeTime output dir
#     outliers.tsv                                 — clock outliers for 03_maple_r2.R
#
# FLAGS:
#   --gene X   run a single gene segment
#   --test     use subsampled data from data/test/

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
  library(dplyr)
})

source(here("revisions", "config.R"))

args      <- commandArgs(trailingOnly = TRUE)
ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
genes     <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}
use_test <- "--test" %in% args

if (use_test) message("*** Running in TEST mode ***")

paths_pruned <- if (use_test) PATHS$test_trees_pruned else PATHS$trees_pruned
paths_meta   <- if (use_test) PATHS$test_metadata     else PATHS$metadata
paths_timed  <- if (use_test) PATHS$test_trees_timed  else PATHS$trees_timed

for (gene in genes) {
  message("\n[", gene, "] TreeTime r1 — clock outlier detection...")

  tree_in  <- file.path(paths_pruned, paste0(gene, "_pruned_r1.nwk"))
  aln_in   <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  meta_in  <- file.path(paths_meta,   paste0(gene, "_metadata.csv"))
  tree_out <- file.path(paths_timed,  paste0(gene, "_timetree_r1.nwk"))
  tt_dir   <- file.path(paths_timed,  paste0(gene, "_treetime_r1_out"))

  dir.create(paths_timed, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(tree_out)) {
    message("  R1 timetree already exists, skipping: ", tree_out)
    next
  }

  for (f in c(tree_in, aln_in, meta_in)) {
    if (!file.exists(f)) { warning("File not found: ", f); next }
  }

  dir.create(tt_dir, recursive = TRUE, showWarnings = FALSE)

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

  tt_files <- list.files(tt_dir,
                         pattern    = "timetree\\.(nwk|nexus|tre)$",
                         full.names = TRUE)
  if (length(tt_files) == 0) {
    warning("TreeTime produced no timetree for ", gene, ". Check: ", tt_dir)
    next
  }

  tt_file  <- tt_files[1]
  timetree <- if (grepl("\\.nexus$", tt_file)) {
    read.nexus(tt_file)
  } else {
    read.tree(tt_file)
  }

  write.tree(timetree, tree_out)
  message("  Written: ", tree_out, " (", Ntip(timetree), " tips)")

  outlier_file <- file.path(tt_dir, "outliers.tsv")
  if (file.exists(outlier_file)) {
    n_out <- nrow(read.delim(outlier_file, stringsAsFactors = FALSE))
    message("  Clock outliers flagged: ", n_out,
            "  (see ", outlier_file, ")")
  } else {
    message("  No outliers.tsv found — TreeTime may not have flagged any.")
  }
}

message("\nDone.")
