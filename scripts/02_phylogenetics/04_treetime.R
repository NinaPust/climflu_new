# scripts/02_phylogenetics/04_treetime.R
#
# PURPOSE: Date the MAPLE tree using TreeTime molecular clock calibration.
#
# INPUT:
#   data/trees/{gene}_tree.nwk
#   data/metadata/{gene}_metadata.csv  — collection_date and accession columns
#
# OUTPUT:
#   data/trees/{gene}_timetree.nwk  — time-calibrated tree (Newick)
#
# REQUIRES: treetime installed and accessible at TOOLS$treetime (see config.R)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
  library(dplyr)
})

source(here("config.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

for (gene in genes) {
  message("\n[", gene, "] TreeTime molecular clock dating...")

  tree_in   <- file.path(PATHS$trees,    paste0(gene, "_tree.nwk"))
  meta_in   <- file.path(PATHS$metadata, paste0(gene, "_metadata.csv"))
  tree_out  <- file.path(PATHS$trees,    paste0(gene, "_timetree.nwk"))
  tt_dir    <- file.path(PATHS$trees,    paste0(gene, "_treetime_out"))

  if (!file.exists(tree_in))  { warning("Tree not found: ", tree_in);  next }
  if (!file.exists(meta_in))  { warning("Metadata not found: ", meta_in); next }

  dir.create(tt_dir, showWarnings = FALSE)

  # Write date file (TSV: name, date) expected by treetime
  meta  <- read_csv(meta_in, show_col_types = FALSE)
  dates <- meta |>
    select(name = accession, date = collection_date) |>
    mutate(date = format(as.Date(date), "%Y-%m-%d"))

  date_file <- file.path(tt_dir, "dates.tsv")
  write_tsv(dates, date_file)

  # Run treetime
  cmd_args <- c(
    "--tree",   tree_in,
    "--dates",  date_file,
    "--outdir", tt_dir,
    "--clock-rate", "auto",
    "--reroot", "least-squares"
  )

  ret <- system2(TOOLS$treetime, args = cmd_args, stdout = TRUE, stderr = TRUE)

  # TreeTime writes timetree.nexus or timetree.nwk depending on version
  tt_files <- list.files(tt_dir, pattern = "timetree\\.(nwk|nexus|tre)$",
                         full.names = TRUE)
  if (length(tt_files) == 0) {
    warning("TreeTime did not produce a timetree output for ", gene,
            ". Check logs in: ", tt_dir)
    next
  }

  timetree <- read.tree(tt_files[1])
  write.tree(timetree, tree_out)
  message("  Written: ", tree_out, " (", Ntip(timetree), " tips)")
}

message("\nDone.")
