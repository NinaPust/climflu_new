# scripts/02_phylogenetics/04_codeml.R
#
# PURPOSE: Estimate pairwise dN and dS for all post-pruning sequences using
#          PAML codeml (runmode=-2). Feeds directly from the MAFFT alignment.
#
# INPUT:
#   data/alignments/{gene}_aligned.fasta         — MAFFT nucleotide alignment
#   data/trees/pruned/{gene}_accessions.txt      — post-MAPLE-r2 accession list
#
# OUTPUT (per batch-pair):
#   data/codeml/{gene}/{pair}/seq.phy            — PHYLIP alignment (alias IDs)
#   data/codeml/{gene}/{pair}/seq_map.csv        — alias → full accession mapping
#   data/codeml/{gene}/{pair}/2ML.dN             — raw codeml dN matrix
#   data/codeml/{gene}/{pair}/2ML.dS             — raw codeml dS matrix
#   data/codeml/{gene}/{pair}/dN_remapped.rds    — dN matrix with full IDs
#   data/codeml/{gene}/{pair}/dS_remapped.rds    — dS matrix with full IDs
#
# DESIGN:
#   Sequences are split into numeric batches (B001, B002, ...).
#   All N*(N+1)/2 pair combinations are run in parallel via mclapply.
#   Short aliases (S00001...) are used in PHYLIP to avoid the 10-char collision
#   problem with EPI_ISL_* accession IDs; names are remapped after codeml.
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(parallel)
})

source(here("config.R"))
source(here("R", "matrix_utils.R"))

args    <- commandArgs(trailingOnly = TRUE)
genes   <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else ALL_GENES
t_start <- proc.time()

# ── Helpers ───────────────────────────────────────────────────────────────────

# Write PHYLIP sequential format with alias IDs (safe for codeml's 10-char limit)
write_phylip <- function(aln_mat, alias_ids, path) {
  char_mat <- as.character(aln_mat)
  n_seq    <- nrow(char_mat)
  n_col    <- ncol(char_mat)
  lines    <- paste(" ", n_seq, n_col)
  for (i in seq_len(n_seq)) {
    nm  <- formatC(alias_ids[i], width = 10, flag = "-")
    seq_str <- paste(char_mat[i, ], collapse = "")
    lines <- c(lines, paste0(nm, "  ", seq_str))
  }
  writeLines(lines, path)
}

make_star_tree <- function(tips) {
  n <- length(tips)
  tree <- list(
    edge        = cbind(rep(n + 1L, n), seq_len(n)),
    tip.label   = tips,
    edge.length = rep(0.1, n),
    Nnode       = 1L
  )
  class(tree) <- "phylo"
  tree
}

codeml_ctl <- function(seqfile, outfile, treefile, model) {
  sprintf(
    "      seqfile  = %s
      outfile  = %s
      treefile = %s

      noisy    = 0
      verbose  = 0
      runmode  = -2

      seqtype  = 1
      CodonFreq = 2
      model    = %s
      NSsites  = 0

      icode    = 0
      fix_kappa = 0
      kappa    = 2
      fix_omega = 0
      omega    = 0.5
      fix_alpha = 1
      alpha    = 0
      Malpha   = 0
      ncatG    = 5
      clock    = 0
      getSE    = 0
      RateAncestor = 0
      cleandata = 1
",
    seqfile, outfile, treefile, model
  )
}

# Run one batch-pair. Returns pair_label on success, NULL on skip, or the label
# with an error attribute on failure.
run_pair <- function(pair, aln_sub, codeml_root, codeml_bin, model, nice) {
  b1 <- pair[1]; b2 <- pair[2]
  pair_label <- paste0(b1, "_", b2)
  pair_dir   <- file.path(codeml_root, pair_label)

  done_file <- file.path(pair_dir, "dS_remapped.rds")
  if (file.exists(done_file)) return(invisible(NULL))  # already done
  t_pair_start <- proc.time()

  dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

  accs <- if (b1 == b2) batches[[b1]] else unique(c(batches[[b1]], batches[[b2]]))
  present <- intersect(accs, rownames(as.matrix(aln_sub)))
  if (length(present) < 4) {
    message("  [", pair_label, "] Too few sequences (", length(present), "), skipping.")
    return(pair_label)
  }

  pair_aln <- as.matrix(aln_sub[present])

  # Assign short aliases to dodge the 10-char PHYLIP limit
  aliases  <- sprintf("S%05d", seq_along(present))
  name_map <- setNames(present, aliases)  # alias -> full accession

  map_df <- data.frame(alias = aliases, accession = present,
                       stringsAsFactors = FALSE)
  write.csv(map_df, file.path(pair_dir, "seq_map.csv"), row.names = FALSE)

  phy_path  <- file.path(pair_dir, "seq.phy")
  tree_path <- file.path(pair_dir, "star.nwk")
  ctl_path  <- file.path(pair_dir, "codeml.ctl")
  out_path  <- file.path(pair_dir, "codeml_out.txt")

  tryCatch({
    write_phylip(pair_aln, aliases, phy_path)
    write.tree(make_star_tree(aliases), tree_path)
    writeLines(codeml_ctl(phy_path, out_path, tree_path, model), ctl_path)

    old_wd <- setwd(pair_dir)
    system2("nice", args = c("-n", nice, codeml_bin, "codeml.ctl"),
            stdout = file.path(pair_dir, "codeml_stdout.txt"),
            stderr = file.path(pair_dir, "codeml_stdout.txt"),
            stdin  = "/dev/null")
    setwd(old_wd)

    ds_raw <- file.path(pair_dir, "2ML.dS")
    dn_raw <- file.path(pair_dir, "2ML.dN")

    if (!file.exists(ds_raw)) {
      message("  [", pair_label, "] codeml produced no output.")
      return(structure(pair_label, error = TRUE))
    }

    saveRDS(remap_paml_matrix(dn_raw, name_map), file.path(pair_dir, "dN_remapped.rds"))
    saveRDS(remap_paml_matrix(ds_raw, name_map), file.path(pair_dir, "dS_remapped.rds"))

    elapsed <- round((proc.time() - t_pair_start)[["elapsed"]])
    message("  [", pair_label, "] Done (", length(present), " seqs, ",
            elapsed, "s).")
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

  # Get post-pruning accession list
  acc_file <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
  if (file.exists(acc_file)) {
    accessions <- readLines(acc_file)
    accessions <- accessions[nchar(trimws(accessions)) > 0]
    message("Loaded ", length(accessions), " post-pruning accessions from ",
            basename(acc_file))
  } else {
    warning("Accession file not found: ", acc_file,
            " — using all sequences in FASTA (test/fallback mode).")
    all_aln    <- ape::read.FASTA(fasta_path)
    accessions <- names(all_aln)
  }

  # Subset alignment to pruned accessions
  if (!exists("all_aln")) all_aln <- ape::read.FASTA(fasta_path)
  aln_sub <- all_aln[intersect(accessions, names(all_aln))]
  rm(all_aln)

  n_seqs <- length(aln_sub)
  message("Sequences after subset: ", n_seqs)
  if (n_seqs < 4) {
    warning("Too few sequences for ", gene, " — skipping.")
    next
  }

  # Build numeric batches: B001, B002, ...
  batch_size <- PARAMS$codeml_batch_size
  n_batches  <- ceiling(n_seqs / batch_size)
  batch_ids  <- sprintf("B%03d", seq_len(n_batches))

  # batches must be visible inside run_pair (via environment capture)
  batches <- setNames(
    lapply(seq_len(n_batches), function(i) {
      names(aln_sub)[((i - 1) * batch_size + 1):min(i * batch_size, n_seqs)]
    }),
    batch_ids
  )

  # All self-pairs + cross-pairs
  cross <- combn(batch_ids, 2, simplify = FALSE)
  self  <- lapply(batch_ids, function(x) c(x, x))
  pairs <- c(self, cross)
  message(n_batches, " batches, ", length(pairs), " pair-jobs.")

  codeml_root <- file.path(PATHS$codeml, gene)
  dir.create(codeml_root, showWarnings = FALSE, recursive = TRUE)

  # Parallel execution
  results <- mclapply(
    pairs,
    function(pair) run_pair(pair, aln_sub, codeml_root,
                             TOOLS$codeml, PARAMS$codeml_model,
                             PARAMS$codeml_nice),
    mc.cores      = PARAMS$codeml_cores,
    mc.preschedule = FALSE   # let workers pick up jobs dynamically
  )

  failed <- Filter(Negate(is.null), results)
  failed <- failed[sapply(failed, function(x) isTRUE(attr(x, "error")))]

  if (length(failed) > 0) {
    fail_log <- file.path(codeml_root, "failed_pairs.txt")
    writeLines(unlist(failed), fail_log)
    message(length(failed), " failed pairs logged to: ", fail_log)
  }

  n_done  <- sum(sapply(results, is.null))
  elapsed <- round((proc.time() - t_gene_start)[["elapsed"]])
  message(gene, ": ", n_done, " pairs completed, ", length(failed),
          " failed. Wall time: ", elapsed, "s (",
          round(elapsed / 60, 1), " min).")
}

total <- round((proc.time() - t_start)[["elapsed"]])
message("\nDone. Total wall time: ", total, "s (", round(total / 60, 1), " min).")
