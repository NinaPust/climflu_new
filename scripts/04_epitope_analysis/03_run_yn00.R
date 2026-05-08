# scripts/04_epitope_analysis/03_run_yn00.R
#
# PURPOSE: Run pairwise yn00 on epitope or non-epitope sub-alignments.
#          Adapted from 02_phylogenetics/04_yn00.R with three key changes:
#            1. batch_size = 500 (fewer pair-jobs, ~35 GB peak vs 385 GB)
#            2. Raw yn00 output files deleted after RDS saved (keeps only
#               dN_remapped.rds / dS_remapped.rds per pair)
#            3. Input/output paths under data/epitope/
#
# INPUT:   data/epitope/{gene}_{type}_aligned.fasta
# OUTPUT:  data/epitope/yn00/{gene}_{type}/{pair}/dN_remapped.rds
#          data/epitope/yn00/{gene}_{type}/failed_pairs.txt
#
# FLAGS:
#   --gene X          run a single gene (default: H1 H3)
#   --type epitope|nonepitope  which sub-alignment (default: both)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(parallel)
})

source(here("config.R"))
source(here("R", "matrix_utils.R"))


args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  c("H1", "H3")
}
types <- if ("--type" %in% args) {
  args[which(args == "--type") + 1]
} else {
  c("epitope", "nonepitope", "ridge", "nonridge")
}

# Larger batches than main pipeline to limit number of pair-jobs and
# peak disk usage (~35 GB vs ~385 GB at batch_size=200).
batch_size <- 500L
epi_dir    <- here("data", "epitope")
t_start    <- proc.time()

# ── Helpers (identical to 04_yn00.R) ─────────────────────────────────────────

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
    cols   <- ((k - 1L) * 3L + 1L):(k * 3L)
    codons <- apply(cm[, cols, drop = FALSE], 1, paste, collapse = "")
    if (!any(codons %in% stop_codons) && !any(cm[, cols] %in% ambig))
      keep_cols[cols] <- TRUE
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
    "      seqfile  = %s\n      outfile  = %s\n      verbose  = 0\n\
      icode    = 0\n      weighting = 0\n      commonf3x4 = 0\n",
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
  t0 <- proc.time()
  dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

  accs    <- if (b1 == b2) batches[[b1]] else
    unique(c(batches[[b1]], batches[[b2]]))
  present <- intersect(accs, rownames(as.matrix(aln_sub)))
  if (length(present) < 4) {
    message("  [", pair_label, "] Too few sequences, skipping.")
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
      data.frame(alias = aliases, accession = present,
                 stringsAsFactors = FALSE),
      file.path(pair_dir, "seq_map.csv"), row.names = FALSE
    )
    writeLines(yn00_ctl(phy_path, out_path), ctl_path)

    old_wd <- setwd(pair_dir)
    system2("nice",
            args   = c("-n", nice, yn00_bin, "yn00.ctl"),
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

    saveRDS(remap_paml_matrix(dn_raw, name_map),
            file.path(pair_dir, "dN_remapped.rds"))
    saveRDS(remap_paml_matrix(ds_raw, name_map),
            file.path(pair_dir, "dS_remapped.rds"))

    # Delete raw yn00 output immediately to cap disk usage.
    # Only the small remapped RDS files are kept until matrix assembly.
    unlink(c(phy_path, ctl_path, out_path,
             file.path(pair_dir, "yn00_stdout.txt"),
             file.path(pair_dir, "seq_map.csv"),
             ds_raw, dn_raw,
             file.path(pair_dir, "2YN.dN"),
             Sys.glob(file.path(pair_dir, "2YN.*")),
             Sys.glob(file.path(pair_dir, "rst*"))))

    elapsed <- round((proc.time() - t0)[["elapsed"]])
    message("  [", pair_label, "] Done (",
            length(present), " seqs, ", elapsed, "s)")
    invisible(NULL)
  }, error = function(e) {
    message("  [", pair_label, "] Error: ", conditionMessage(e))
    structure(pair_label, error = TRUE)
  })
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for (gene in genes) {
  for (type in types) {
    tag <- paste0(gene, "_", type)
    message("\n========== ", tag, " ==========")

    fasta_path <- file.path(epi_dir, paste0(tag, "_aligned.fasta"))
    if (!file.exists(fasta_path)) {
      warning("FASTA not found, skipping: ", fasta_path); next
    }

    aln_sub  <- read.FASTA(fasta_path)
    n_seqs   <- length(aln_sub)
    message("Sequences: ", n_seqs)
    if (n_seqs < 4) {
      warning("Too few sequences for ", tag, " — skipping."); next
    }

    n_batches <- ceiling(n_seqs / batch_size)
    batch_ids <- sprintf("B%03d", seq_len(n_batches))
    batches   <- setNames(
      lapply(seq_len(n_batches), function(i) {
        names(aln_sub)[((i - 1) * batch_size + 1):
                         min(i * batch_size, n_seqs)]
      }),
      batch_ids
    )

    pairs <- c(
      lapply(batch_ids, function(x) c(x, x)),
      combn(batch_ids, 2, simplify = FALSE)
    )
    message(n_batches, " batches, ", length(pairs), " pair-jobs")

    yn00_root <- file.path(epi_dir, "yn00", tag)
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
      message(length(failed), " failed pairs logged: ", fail_log)
    }

    elapsed <- round((proc.time() - t_start)[["elapsed"]])
    message(tag, ": ", sum(sapply(results, is.null)), " pairs done, ",
            length(failed), " failed. Wall time: ",
            round(elapsed / 60, 1), " min")
  }
}

message("\nDone.")
