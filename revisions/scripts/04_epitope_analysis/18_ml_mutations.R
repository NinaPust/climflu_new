# revisions/scripts/04_epitope_analysis/18_ml_mutations.R
#
# PURPOSE: Identify individual non-synonymous mutation events on TERMINAL
#          branches using TreeTime's maximum-likelihood ancestral sequences,
#          replacing the parsimony reconstruction in 09_substitutions.R as the
#          source for figure 2.
#
# WHY ML: 09_substitutions.R uses phangorn::ancestral.pars with MPR. Where a
#   site has several equally parsimonious ancestral states, that code resolves
#   the tie with `which(x > 0)[1]` — i.e. whichever state happens to come first
#   in the level ordering. ML assigns the ancestral state under a fitted
#   substitution model instead, so ties are resolved by evidence rather than by
#   ordering. On short terminal branches with dense sampling the two will often
#   agree; this script reports the overlap against the parsimony output so the
#   size of the difference is a measured quantity rather than an assumption.
#
# WHY NOT SLAC: SLAC's codon-model reconstruction would in principle be better
#   still, but its JSON publishes internal-node states under names (Node1..NodeN)
#   that cannot be mapped back to tree positions — the tree it exports carries
#   tip labels only, and it collapses zero-length internal branches (3898 -> 3276
#   for H1), so the input tree does not map either. TreeTime's timetree DOES
#   carry internal node labels and its ancestral FASTA is keyed by them, so the
#   tip -> parent lookup here is exact.
#
# TERMINAL BRANCHES ONLY. A mutation is recorded against the tip it occurred on,
#   because figure 2 tests each mutation against that sample's climate, and
#   climate is a property of a sample rather than of an internal node.
#
# OUTPUT SCHEMA matches 09_substitutions.R exactly, so 19_wilcoxon_ml.R (and
#   scripts/05_final_figs/fig2_mutations.R) consume it unchanged:
#     accession, gene, region, residue_pos, parent_aa, tip_aa,
#     parent_codon, tip_codon, change
#
# INPUT:
#   data/trees/timed/{gene}_timetree.nwk                   internal node labels
#   data/trees/timed/{gene}_treetime_out/ancestral_sequences.fasta
#   revisions/data/epitope/{gene}_site_map.csv             codon -> region
#
# OUTPUT:
#   revisions/data/epitope/{gene}_{region}_mutations_ml.csv
#   revisions/data/epitope/{gene}_ml_mutation_summary.csv   per-region totals and
#                                                           overlap with parsimony
#
# FLAGS:
#   --gene X       run a single gene (default: H1 H3)
#   --region X     ridge | epitope | HA1 | HA2 | full (default: ridge epitope)
#   --no-compare   skip the comparison against the parsimony mutations
#
# RUNTIME: a few minutes per gene; the work is 565 vectorised codon comparisons
#   across ~11k tips. Memory peak is the alignment matrix, roughly 250 MB.

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(Biostrings)
  library(dplyr)
  library(readr)
  library(tibble)
})

source(here("revisions", "config.R"))

args       <- commandArgs(trailingOnly = TRUE)
genes      <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
regions    <- if ("--region" %in% args) args[which(args == "--region") + 1] else
              c("ridge", "epitope")
do_compare <- !("--no-compare" %in% args)

GENETIC <- Biostrings::GENETIC_CODE
epi_dir <- PATHS$epitope
dir.create(epi_dir, recursive = TRUE, showWarnings = FALSE)


#' Map each CDS codon to its region and mature residue position.
#'
#' The CDS columns are the site map's alignment columns in sorted order, so
#' codon k occupies positions (3k-2):(3k) of that vector. Region and residue
#' position are read from the codon's first column; all three columns of a codon
#' share both by construction.
codon_table <- function(site_map_path) {
  sm   <- read_csv(site_map_path, show_col_types = FALSE)
  cols <- sort(unique(sm$alignment_col))
  if (length(cols) %% 3 != 0)
    stop("site map column count (", length(cols), ") is not a multiple of 3")
  first <- cols[seq(1, length(cols), by = 3)]
  idx   <- match(first, sm$alignment_col)
  tibble(
    codon       = seq_along(first),
    region_raw  = sm$region[idx],
    residue_pos = sm$residue_pos[idx]
  ) |>
    mutate(
      is_ridge   = grepl("_ridge$", region_raw),
      is_epitope = grepl("^epitope", region_raw),
      is_HA1     = region_raw == "HA1",
      is_HA2     = grepl("^HA2", region_raw)
    ) -> tab
  attr(tab, "cols") <- cols
  tab
}

#' Codon indices belonging to a named region.
region_codons <- function(tab, region) {
  switch(region,
    full    = tab$codon,
    ridge   = tab$codon[tab$is_ridge],
    epitope = tab$codon[tab$is_epitope],
    HA1     = tab$codon[tab$is_HA1],
    HA2     = tab$codon[tab$is_HA2],
    stop("unknown region: ", region)
  )
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for (gene in genes) {
  message("\n========== ", gene, " ==========")

  tree_path <- file.path(here(), "data", "trees", "timed",
                         paste0(gene, "_timetree.nwk"))
  anc_path  <- file.path(here(), "data", "trees", "timed",
                         paste0(gene, "_treetime_out"), "ancestral_sequences.fasta")
  smap_path <- file.path(epi_dir, paste0(gene, "_site_map.csv"))
  miss <- c(tree_path, anc_path, smap_path)
  miss <- miss[!file.exists(miss)]
  if (length(miss)) { warning(gene, ": missing ", paste(basename(miss), collapse = ", ")); next }

  tr <- read.tree(tree_path)
  if (is.null(tr$node.label))
    stop(gene, ": timetree has no internal node labels — cannot identify parents")

  tab  <- codon_table(smap_path)
  cols <- attr(tab, "cols")
  message("  CDS: ", length(cols), " nt = ", nrow(tab), " codons")

  anc <- readDNAStringSet(anc_path)
  if (unique(width(anc))[1] < max(cols))
    stop(gene, ": ancestral sequences are narrower than the site map requires")

  # Tip -> parent label, from the labelled timetree. Exact, no inference.
  tip_edge <- match(seq_len(Ntip(tr)), tr$edge[, 2])
  par_idx  <- tr$edge[tip_edge, 1]
  par_lab  <- tr$node.label[par_idx - Ntip(tr)]
  tips     <- tr$tip.label

  have <- tips %in% names(anc) & par_lab %in% names(anc)
  if (!all(have))
    message("  ", sum(!have), " of ", length(tips),
            " tips dropped (tip or parent sequence absent from ancestral FASTA)")
  tips <- tips[have]; par_lab <- par_lab[have]
  message("  terminal branches to scan: ", length(tips))

  A  <- as.matrix(anc)
  TM <- A[tips,    cols, drop = FALSE]   # tip sequences, CDS columns
  PM <- A[par_lab, cols, drop = FALSE]   # their parents' ML sequences
  rm(A); invisible(gc(verbose = FALSE))

  # ── Scan codons ─────────────────────────────────────────────────────────────
  # One vectorised pass per codon over all terminal branches.
  t0  <- proc.time()
  hits <- vector("list", nrow(tab))
  for (k in seq_len(nrow(tab))) {
    j  <- (k - 1L) * 3L
    tc <- toupper(paste0(TM[, j + 1L], TM[, j + 2L], TM[, j + 3L]))
    pc <- toupper(paste0(PM[, j + 1L], PM[, j + 2L], PM[, j + 3L]))
    d  <- which(tc != pc)
    if (!length(d)) next

    ta <- GENETIC[tc[d]]; pa <- GENETIC[pc[d]]
    # Drop ambiguous or gapped codons (NA) and any stop codon: neither is a
    # scoreable amino-acid change.
    ok <- !is.na(ta) & !is.na(pa) & ta != "*" & pa != "*" & ta != pa
    if (!any(ok)) next

    hits[[k]] <- tibble(
      accession    = tips[d][ok],
      codon_idx    = k,
      residue_pos  = tab$residue_pos[k],
      parent_codon = pc[d][ok],
      tip_codon    = tc[d][ok],
      parent_aa    = unname(pa[ok]),
      tip_aa       = unname(ta[ok])
    )
  }
  muts <- bind_rows(hits)
  message("  non-synonymous events on terminal branches: ", nrow(muts),
          "   [", round((proc.time() - t0)[["elapsed"]]), "s]")

  # ── Write per region ────────────────────────────────────────────────────────
  summary_rows <- list()
  for (r in regions) {
    rc <- region_codons(tab, r)
    mr <- muts |>
      filter(codon_idx %in% rc) |>
      mutate(gene = gene, region = r,
             change = paste0(parent_aa, residue_pos, tip_aa)) |>
      select(accession, gene, region, residue_pos, parent_aa, tip_aa,
             parent_codon, tip_codon, change) |>
      arrange(residue_pos, change)

    out_path <- file.path(epi_dir, paste0(gene, "_", r, "_mutations_ml.csv"))
    write_csv(mr, out_path)
    message("  [", r, "] ", nrow(mr), " events, ",
            n_distinct(mr$change), " distinct changes -> ", basename(out_path))

    row <- tibble(gene = gene, region = r, codons = length(rc),
                  n_events_ml = nrow(mr), n_changes_ml = n_distinct(mr$change))

    # ── Compare against the parsimony output, where it exists ────────────────
    if (do_compare) {
      old <- file.path(here(), "data", "epitope",
                       paste0(gene, "_", r, "_mutations.csv"))
      if (file.exists(old)) {
        op <- read_csv(old, show_col_types = FALSE)
        key_ml  <- paste(mr$accession, mr$change)
        key_par <- paste(op$accession, op$change)
        row$n_events_parsimony <- nrow(op)
        row$shared   <- length(intersect(key_ml, key_par))
        row$ml_only  <- length(setdiff(key_ml, key_par))
        row$par_only <- length(setdiff(key_par, key_ml))
        row$pct_agreement <- round(100 * row$shared /
                                   length(union(key_ml, key_par)), 1)
        message("      vs parsimony: shared ", row$shared,
                "  ML-only ", row$ml_only, "  parsimony-only ", row$par_only,
                "  (", row$pct_agreement, "% agreement)")
      } else {
        message("      (no parsimony file at ", basename(old), " to compare)")
      }
    }
    summary_rows[[r]] <- row
  }

  sum_path <- file.path(epi_dir, paste0(gene, "_ml_mutation_summary.csv"))
  write_csv(bind_rows(summary_rows), sum_path)
  message("  Written: ", basename(sum_path))

  # Most recurrent changes, as a sanity check on the output.
  for (r in regions) {
    f <- file.path(epi_dir, paste0(gene, "_", r, "_mutations_ml.csv"))
    if (!file.exists(f)) next
    top <- read_csv(f, show_col_types = FALSE) |>
      count(change, sort = TRUE) |> head(6)
    if (nrow(top))
      message("  [", r, "] most recurrent: ",
              paste(sprintf("%s (%d)", top$change, top$n), collapse = ", "))
  }
}

message("\nDone.")
