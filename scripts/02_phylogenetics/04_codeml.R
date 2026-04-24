# scripts/02_phylogenetics/02_codeml.R
#
# PURPOSE: Run PAML codeml on each pairwise MACSE-aligned bin to estimate
#          pairwise dN and dS values.
#
# INPUT:
#   data/alignments/{gene}/macse/{pair}_aligned.rds  (from 02_macse.R)
#
# OUTPUT:
#   data/codeml/{gene}/{pair}/2ML.dN  — pairwise dN matrix
#   data/codeml/{gene}/{pair}/2ML.dS  — pairwise dS matrix

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("config.R"))
source(here("R", "alignment_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else GENES

# codeml control file template
codeml_ctl_template <- '
      seqfile = {seqfile}
      outfile = {outfile}
      treefile = {treefile}

      noisy = 0
      verbose = 0
      runmode = -2    * pairwise distances

      seqtype = 1     * codons
      CodonFreq = 2   * F3x4
      model = {model}
      NSsites = 0

      icode = 0
      fix_kappa = 0
      kappa = 2
      fix_omega = 0
      omega = 0.5
      fix_alpha = 1
      alpha = 0
      Malpha = 0
      ncatG = 5
      clock = 0
      getSE = 0
      RateAncestor = 0
      cleandata = 1
'

write_ctl <- function(seqfile, outfile, treefile, model, ctl_path) {
  ctl <- codeml_ctl_template
  ctl <- gsub("\\{seqfile\\}",  seqfile,  ctl)
  ctl <- gsub("\\{outfile\\}",  outfile,  ctl)
  ctl <- gsub("\\{treefile\\}", treefile, ctl)
  ctl <- gsub("\\{model\\}",    model,    ctl)
  writeLines(ctl, ctl_path)
}

for (gene in genes) {
  message("\n[", gene, "] codeml dN/dS...")

  macse_dir   <- file.path(PATHS$alignments, gene, "macse")
  codeml_root <- file.path(PATHS$codeml,    gene)
  dir.create(codeml_root, showWarnings = FALSE, recursive = TRUE)

  rds_files <- list.files(macse_dir, pattern = "_aligned\\.rds$", full.names = TRUE)
  if (length(rds_files) == 0) {
    warning("No MACSE RDS files found for gene ", gene, ". Run 01_alignment/02_macse.R first.")
    next
  }

  failed <- character(0)

  for (rds_file in rds_files) {
    pair_label <- sub("_aligned\\.rds$", "", basename(rds_file))
    pair_dir   <- file.path(codeml_root, pair_label)

    if (file.exists(file.path(pair_dir, "2ML.dS"))) {
      message("  [", pair_label, "] Already done, skipping.")
      next
    }

    dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

    bin <- tryCatch(readRDS(rds_file), error = function(e) NULL)
    if (is.null(bin)) {
      message("  [", pair_label, "] Could not read RDS, skipping.")
      failed <- c(failed, pair_label)
      next
    }

    # Write alignment as PHYLIP (required by codeml)
    phy_path <- file.path(pair_dir, "seq.phy")
    tryCatch(
      write.phy(bin, phy_path),
      error = function(e) {
        # ape::write.phy may not exist in all versions; fall back to FASTA
        write.dna(bin, phy_path, format = "interleaved")
      }
    )

    # Minimal star tree (codeml runmode=-2 doesn't use topology)
    n_taxa <- nrow(as.matrix(bin))
    tip_labels <- rownames(as.matrix(bin))
    star_tree <- strtree(tip_labels)
    tree_path <- file.path(pair_dir, "star.nwk")
    write.tree(star_tree, tree_path)

    ctl_path <- file.path(pair_dir, "codeml.ctl")
    out_path <- file.path(pair_dir, "codeml_out.txt")
    write_ctl(phy_path, out_path, tree_path, PARAMS$codeml_model, ctl_path)

    # Run codeml from pair_dir (codeml writes output relative to cwd)
    old_wd <- setwd(pair_dir)
    ret <- system2(TOOLS$codeml, args = "codeml.ctl", stdout = TRUE, stderr = TRUE)
    setwd(old_wd)

    if (!file.exists(file.path(pair_dir, "2ML.dS"))) {
      message("  [", pair_label, "] codeml did not produce output.")
      failed <- c(failed, pair_label)
    } else {
      message("  [", pair_label, "] Done.")
    }
  }

  if (length(failed) > 0) {
    fail_log <- file.path(codeml_root, "failed_pairs.txt")
    writeLines(failed, fail_log)
    message("  ", length(failed), " failed pairs logged: ", fail_log)
  }
}

message("\nDone.")


# ── Helper: make a star tree from tip labels ──────────────────────────────────
strtree <- function(tips) {
  n <- length(tips)
  tree <- list(
    edge       = cbind(rep(n + 1L, n), seq_len(n)),
    tip.label  = tips,
    edge.length = rep(0.1, n),
    Nnode      = 1L
  )
  class(tree) <- "phylo"
  tree
}
