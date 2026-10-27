# scripts/02_phylogenetics/04z_test_yn00.R
#
# PURPOSE: Smoke-test yn00 on a 120-sequence subset of H1 to verify that
#          the binary works, stop-codon stripping is correct, alias mapping
#          round-trips, and 2YN.dS output is produced and parseable.
#
# OUTPUT:  data/test/yn00_smoke/H1/{pair}/

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(parallel)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "matrix_utils.R"))

GENE       <- "H1"
N_SUBSET   <- 120
BATCH_SIZE <- 50L
N_CORES    <- 3L
OUT_ROOT   <- here("revisions", "data", "test", "yn00_smoke", GENE)
FASTA_IN   <- file.path(PATHS$alignments, paste0(GENE, "_aligned.fasta"))
YN00_BIN   <- "/home/nina/paml4.8/bin/yn00"

t0 <- proc.time()

# ── 1. Load + subsample ───────────────────────────────────────────────────────
message("Reading alignment: ", FASTA_IN)
aln_all <- ape::read.FASTA(FASTA_IN)
message("Total sequences: ", length(aln_all))

set.seed(42)
aln_sub <- aln_all[sample(length(aln_all), min(N_SUBSET, length(aln_all)))]
message("Subsampled to: ", length(aln_sub))

# ── 2. Batches + pairs ────────────────────────────────────────────────────────
n_seqs    <- length(aln_sub)
n_batches <- ceiling(n_seqs / BATCH_SIZE)
batch_ids <- sprintf("B%03d", seq_len(n_batches))
batches   <- setNames(
  lapply(seq_len(n_batches), function(i)
    names(aln_sub)[((i - 1) * BATCH_SIZE + 1):min(i * BATCH_SIZE, n_seqs)]),
  batch_ids)
pairs <- c(lapply(batch_ids, function(x) c(x, x)),
           combn(batch_ids, 2, simplify = FALSE))
message("Batches: ", n_batches, "  Pairs: ", length(pairs))

dir.create(OUT_ROOT, showWarnings = FALSE, recursive = TRUE)

# ── Helpers ───────────────────────────────────────────────────────────────────

# Prepare alignment for yn00:
#   1. Auto-detect reading frame offset (0/1/2) that minimises stop codons.
#   2. Remove any codon column where any sequence has a stop codon or gap.
# This replicates codeml's cleandata=1 behaviour, which yn00 lacks.
clean_for_yn00 <- function(cm) {
  stop_codons <- c("taa", "tag", "tga")
  ambig       <- c("-", "n", "?")

  # Pick the offset with the fewest stop-codon columns
  count_stops <- function(offset) {
    cols_start <- seq(offset + 1, ncol(cm) - 2, by = 3)
    sum(sapply(cols_start, function(k) {
      any(apply(cm[, k:(k + 2), drop = FALSE], 1,
                function(x) paste(x, collapse = "") %in% stop_codons)
      )
    }))
  }
  best_offset <- which.min(sapply(0:2, count_stops)) - 1L
  cm <- cm[, (best_offset + 1):ncol(cm), drop = FALSE]

  # Remove codon columns with any stop codon or gap in any sequence
  n_codons   <- floor(ncol(cm) / 3)
  keep_cols  <- logical(ncol(cm))
  for (k in seq_len(n_codons)) {
    cols   <- ((k - 1L) * 3L + 1L):(k * 3L)
    codons <- apply(cm[, cols, drop = FALSE], 1, paste, collapse = "")
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

# ── 3. Run pairs in parallel ──────────────────────────────────────────────────
run_pair <- function(pair) {
  t1    <- proc.time()
  b1    <- pair[1]; b2 <- pair[2]
  label    <- paste0(b1, "_", b2)
  pair_dir <- file.path(OUT_ROOT, label)

  if (file.exists(file.path(pair_dir, "dS_remapped.rds"))) {
    message("  [", label, "] Already done, skipping.")
    return(invisible(NULL))
  }

  dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

  accs    <- if (b1 == b2) batches[[b1]] else unique(c(batches[[b1]], batches[[b2]]))
  present <- intersect(accs, rownames(as.matrix(aln_sub)))
  if (length(present) < 4) {
    message("  [", label, "] Too few sequences, skipping.")
    return(label)
  }

  pm      <- as.matrix(aln_sub[present])
  aliases <- sprintf("S%05d", seq_along(present))
  nm      <- setNames(present, aliases)

  write_phylip(pm, aliases, file.path(pair_dir, "seq.phy"))
  write.csv(
    data.frame(alias = aliases, accession = present, stringsAsFactors = FALSE),
    file.path(pair_dir, "seq_map.csv"), row.names = FALSE)

  ctl <- sprintf(
    "      seqfile  = %s\n      outfile  = %s\n      verbose  = 0\n      icode    = 0\n      weighting = 0\n      commonf3x4 = 0\n",
    file.path(pair_dir, "seq.phy"),
    file.path(pair_dir, "yn00_out.txt"))
  writeLines(ctl, file.path(pair_dir, "yn00.ctl"))

  tryCatch({
    old_wd <- setwd(pair_dir)
    system2("nice", c("-n", "10", YN00_BIN, "yn00.ctl"),
            stdout = "stdout.txt", stderr = "stdout.txt", stdin = "/dev/null")
    setwd(old_wd)

    if (!file.exists(file.path(pair_dir, "2YN.dS"))) {
      message("  [", label, "] yn00 produced no output.")
      return(structure(label, error = TRUE))
    }

    saveRDS(remap_paml_matrix(file.path(pair_dir, "2YN.dN"), nm),
            file.path(pair_dir, "dN_remapped.rds"))
    saveRDS(remap_paml_matrix(file.path(pair_dir, "2YN.dS"), nm),
            file.path(pair_dir, "dS_remapped.rds"))

    elapsed <- round((proc.time() - t1)[["elapsed"]])
    message("  [", label, "] Done (", length(present), " seqs, ", elapsed, "s).")
    invisible(NULL)
  }, error = function(e) {
    message("  [", label, "] Error: ", conditionMessage(e))
    structure(label, error = TRUE)
  })
}

results <- mclapply(pairs, run_pair, mc.cores = N_CORES, mc.preschedule = FALSE)
failed  <- Filter(function(x) isTRUE(attr(x, "error")), results)

# ── 4. Report ─────────────────────────────────────────────────────────────────
pair_dirs <- list.dirs(OUT_ROOT, full.names = TRUE, recursive = FALSE)
pair_dirs <- pair_dirs[file.exists(file.path(pair_dirs, "dS_remapped.rds"))]

message("\n--- Results ---")
message("Completed: ", length(pair_dirs), " / ", length(pairs), " pairs")
message("Failed:    ", length(failed))

if (length(pair_dirs) > 0) {
  m <- readRDS(file.path(pair_dirs[1], "dS_remapped.rds"))
  message("Matrix dim:        ", nrow(m), " x ", ncol(m))
  message("Full accession IDs: ", all(grepl("^EPI_ISL_", rownames(m))))
  message("Name collisions:    ", if (anyDuplicated(rownames(m)) == 0) "NONE" else "YES")
}

total <- round((proc.time() - t0)[["elapsed"]])
message("Total time: ", total, "s")
message("\nDone.")
