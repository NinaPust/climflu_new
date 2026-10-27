# scripts/02_phylogenetics/04_yn00.R
#
# PURPOSE: Estimate pairwise dN and dS for all post-pruning sequences using
#          PAML yn00 (Yang & Nielsen 2000 method). Much faster than codeml
#          runmode=-2 (~100-200x) at the cost of using an approximate rather
#          than full-ML estimator. Suitable for downstream phylolm analysis.
#
# INPUT:
#   data/alignments/{gene}_aligned.fasta         — MAFFT nucleotide alignment
#   data/trees/pruned/{gene}_accessions.txt      — post-MAPLE-r2 accession list
#
# OUTPUT (per batch-pair):
#   data/yn00/{gene}/{pair}/seq.phy              — PHYLIP alignment (alias IDs)
#   data/yn00/{gene}/{pair}/seq_map.csv          — alias -> full accession mapping
#   data/yn00/{gene}/{pair}/2YN.dN               — raw yn00 dN matrix
#   data/yn00/{gene}/{pair}/2YN.dS               — raw yn00 dS matrix
#   data/yn00/{gene}/{pair}/dN_remapped.rds      — dN matrix with full IDs
#   data/yn00/{gene}/{pair}/dS_remapped.rds      — dS matrix with full IDs
#
# DESIGN:
#   Sequences are split into numeric batches (B001, B002, ...).
#   All N*(N+1)/2 pair combinations are run in parallel via mclapply.
#   Short aliases (S00001...) are used in PHYLIP to avoid the 10-char collision
#   problem with EPI_ISL_* accession IDs; names are remapped after yn00.
#   The alignment is cleaned per-pair (frame offset + gap/stop removal) because
#   yn00 lacks codeml's cleandata option.
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(parallel)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "matrix_utils.R"))

args    <- commandArgs(trailingOnly = TRUE)
genes   <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES
t_start <- proc.time()

# ── Helpers ───────────────────────────────────────────────────────────────────

# Clean alignment for yn00:
#   1. Auto-detect reading frame offset (0/1/2) minimising stop-codon columns.
#   2. Remove codon columns with any stop codon or gap in any sequence.
# yn00 has no cleandata option, so we replicate codeml cleandata=1 here.
clean_for_yn00 <- function(cm) {
  stop_codons <- c("taa", "tag", "tga")
  ambig       <- c("-", "n", "?")

  count_stop_cols <- function(offset) {
    starts <- seq(offset + 1, ncol(cm) - 2, by = 3)
    sum(sapply(starts, function(k) {
      any(apply(cm[, k:(k + 2), drop = FALSE], 1,
                function(x) paste(x, collapse = "") %in% stop_codons))
    }))
  }
  best_offset <- which.min(sapply(0:2, count_stop_cols)) - 1L
  cm <- cm[, (best_offset + 1):ncol(cm), drop = FALSE]

  n_codons  <- floor(ncol(cm) / 3)
  keep_cols <- logical(ncol(cm))
  for (k in seq_len(n_codons)) {
    cols <- ((k - 1L) * 3L + 1L):(k * 3L)
    codons   <- apply(cm[, cols, drop = FALSE], 1, paste, collapse = "")
    has_stop <- any(codons %in% stop_codons)
    has_gap  <- any(cm[, cols] %in% ambig)
    if (!has_stop && !has_gap) keep_cols[cols] <- TRUE
  }
  cm[, keep_cols, drop = FALSE]
}

write_phylip <- function(aln_mat, alias_ids, path) {
  cm    <- clean_for_yn00(as.character(aln_mat))
  lines <- paste(" ", nrow(cm), ncol(cm))
  for (i in seq_len(nrow(cm))) {
    nm    <- formatC(alias_ids[i], width = 10, flag = "-")
    lines <- c(lines, paste0(nm, "  ", paste(cm[i, ], collapse = "")))
  }
  writeLines(lines, path)
}

yn00_ctl <- function(seqfile, outfile) {
  sprintf(
    "      seqfile  = %s
      outfile  = %s
      verbose  = 0
      icode    = 0
      weighting = 0
      commonf3x4 = 0
",
    seqfile, outfile
  )
}

run_pair <- function(pair, aln_sub, yn00_root, yn00_bin, nice) {
  b1 <- pair[1]; b2 <- pair[2]
  pair_label <- paste0(b1, "_", b2)
  pair_dir   <- file.path(yn00_root, pair_label)

  if (file.exists(file.path(pair_dir, "dS_remapped.rds"))) {
    return(invisible(NULL))
  }
  t_pair_start <- proc.time()

  dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

  accs    <- if (b1 == b2) batches[[b1]] else unique(c(batches[[b1]], batches[[b2]]))
  present <- intersect(accs, rownames(as.matrix(aln_sub)))
  if (length(present) < 4) {
    message("  [", pair_label, "] Too few sequences (", length(present), "), skipping.")
    return(pair_label)
  }

  pair_aln <- as.matrix(aln_sub[present])
  aliases  <- sprintf("S%05d", seq_along(present))
  name_map <- setNames(present, aliases)

  phy_path <- file.path(pair_dir, "seq.phy")
  ctl_path <- file.path(pair_dir, "yn00.ctl")
  out_path <- file.path(pair_dir, "yn00_out.txt")

  tryCatch({
    write_phylip(pair_aln, aliases, phy_path)
    write.csv(
      data.frame(alias = aliases, accession = present, stringsAsFactors = FALSE),
      file.path(pair_dir, "seq_map.csv"), row.names = FALSE)
    writeLines(yn00_ctl(phy_path, out_path), ctl_path)

    old_wd <- setwd(pair_dir)
    system2("nice", args = c("-n", nice, yn00_bin, "yn00.ctl"),
            stdout = file.path(pair_dir, "yn00_stdout.txt"),
            stderr = file.path(pair_dir, "yn00_stdout.txt"),
            stdin  = "/dev/null")
    setwd(old_wd)

    ds_raw <- file.path(pair_dir, "2YN.dS")
    dn_raw <- file.path(pair_dir, "2YN.dN")

    if (!file.exists(ds_raw)) {
      message("  [", pair_label, "] yn00 produced no output.")
      return(structure(pair_label, error = TRUE))
    }

    saveRDS(remap_paml_matrix(dn_raw, name_map), file.path(pair_dir, "dN_remapped.rds"))
    saveRDS(remap_paml_matrix(ds_raw, name_map), file.path(pair_dir, "dS_remapped.rds"))

    elapsed <- round((proc.time() - t_pair_start)[["elapsed"]])
    message("  [", pair_label, "] Done (", length(present), " seqs, ", elapsed, "s).")
    invisible(NULL)
  }, error = function(e) {
    message("  [", pair_label, "] Error: ", conditionMessage(e))
    structure(pair_label, error = TRUE)
  })
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for (gene in genes) {
  t_gene_start <- proc.time()
  message("\n========== ", gene, " ==========")

  fasta_path <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  if (!file.exists(fasta_path)) {
    warning("FASTA not found: ", fasta_path, " — skipping ", gene)
    next
  }

  acc_file <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
  if (file.exists(acc_file)) {
    accessions <- readLines(acc_file)
    accessions <- accessions[nchar(trimws(accessions)) > 0]
    message("Loaded ", length(accessions), " post-pruning accessions from ",
            basename(acc_file))
  } else {
    warning("Accession file not found: ", acc_file,
            " — using all sequences in FASTA (fallback mode).")
    all_aln    <- ape::read.FASTA(fasta_path)
    accessions <- names(all_aln)
  }

  if (!exists("all_aln")) all_aln <- ape::read.FASTA(fasta_path)
  aln_sub <- all_aln[intersect(accessions, names(all_aln))]
  rm(all_aln)

  n_seqs <- length(aln_sub)
  message("Sequences after subset: ", n_seqs)
  if (n_seqs < 4) {
    warning("Too few sequences for ", gene, " — skipping.")
    next
  }

  batch_size <- PARAMS$codeml_batch_size
  n_batches  <- ceiling(n_seqs / batch_size)
  batch_ids  <- sprintf("B%03d", seq_len(n_batches))

  batches <- setNames(
    lapply(seq_len(n_batches), function(i) {
      names(aln_sub)[((i - 1) * batch_size + 1):min(i * batch_size, n_seqs)]
    }),
    batch_ids
  )

  cross <- combn(batch_ids, 2, simplify = FALSE)
  self  <- lapply(batch_ids, function(x) c(x, x))
  pairs <- c(self, cross)
  message(n_batches, " batches, ", length(pairs), " pair-jobs.")

  yn00_root <- file.path(PATHS$yn00, gene)
  dir.create(yn00_root, showWarnings = FALSE, recursive = TRUE)

  results <- mclapply(
    pairs,
    function(pair) run_pair(pair, aln_sub, yn00_root,
                            TOOLS$yn00, PARAMS$codeml_nice),
    mc.cores       = PARAMS$codeml_cores,
    mc.preschedule = FALSE
  )

  failed <- Filter(function(x) isTRUE(attr(x, "error")), results)

  if (length(failed) > 0) {
    fail_log <- file.path(yn00_root, "failed_pairs.txt")
    writeLines(unlist(failed), fail_log)
    message(length(failed), " failed pairs logged to: ", fail_log)
  }

  n_done  <- sum(sapply(results, is.null))
  elapsed <- round((proc.time() - t_gene_start)[["elapsed"]])
  message(gene, ": ", n_done, " pairs completed, ", length(failed),
          " failed. Wall time: ", elapsed, "s (", round(elapsed / 60, 1), " min).")
}

total <- round((proc.time() - t_start)[["elapsed"]])
message("\nDone. Total wall time: ", total, "s (", round(total / 60, 1), " min).")
