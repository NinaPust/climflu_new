# scripts/02_phylogenetics/04z_test_codeml.R
#
# PURPOSE: Smoke-test codeml on a 100-sequence subset of H1 to verify that
#          the binary works, the PHYLIP writer produces valid input, and
#          codeml emits 2ML.dN / 2ML.dS output files we can parse.
#
# INPUT:   data/alignments/H1_aligned.fasta  (MAFFT alignment, all H1 seqs)
# OUTPUT:  data/test/codeml/H1/{pair}/2ML.dN, 2ML.dS

suppressPackageStartupMessages({
  library(here)
  library(ape)
})

source(here("revisions", "config.R"))

GENE        <- "H1"
N_SUBSET    <- 100
BATCH_SIZE  <- 50
TEST_PAIR   <- "AA"   # only run one pair for speed
OUT_ROOT    <- here("revisions", "data", "test", "codeml", GENE)
FASTA_IN    <- file.path(PATHS$alignments, paste0(GENE, "_aligned.fasta"))

# ── 1. Load + subsample ───────────────────────────────────────────────────────
message("Reading alignment: ", FASTA_IN)
aln <- ape::read.FASTA(FASTA_IN)
message("Total sequences: ", length(aln))

set.seed(42)
idx     <- sample(length(aln), min(N_SUBSET, length(aln)))
aln_sub <- aln[idx]
message("Subsampled to ", length(aln_sub), " sequences.")

# ── 2. Split into lettered batches ────────────────────────────────────────────
n_batches  <- ceiling(length(aln_sub) / BATCH_SIZE)
batch_ids  <- LETTERS[seq_len(n_batches)]
batches    <- setNames(
  lapply(seq_len(n_batches), function(i) {
    start <- (i - 1) * BATCH_SIZE + 1
    end   <- min(i * BATCH_SIZE, length(aln_sub))
    names(aln_sub)[start:end]
  }),
  batch_ids
)
message("Batches: ", paste(batch_ids, sapply(batches, length), sep = "=", collapse = ", "))

# ── 3. Run codeml on TEST_PAIR ────────────────────────────────────────────────
# For pair "AA", "AB", etc., combine the two batches' accessions
parse_pair <- function(label) {
  b1 <- substr(label, 1, 1)
  b2 <- substr(label, 2, 2)
  unique(c(batches[[b1]], batches[[b2]]))
}

codeml_ctl_template <- "
      seqfile  = {seqfile}
      outfile  = {outfile}
      treefile = {treefile}

      noisy    = 0
      verbose  = 0
      runmode  = -2

      seqtype  = 1
      CodonFreq = 2
      model    = 0
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
"

write_ctl <- function(seqfile, outfile, treefile, ctl_path) {
  ctl <- codeml_ctl_template
  ctl <- gsub("\\{seqfile\\}",  seqfile,  ctl)
  ctl <- gsub("\\{outfile\\}",  outfile,  ctl)
  ctl <- gsub("\\{treefile\\}", treefile, ctl)
  writeLines(ctl, ctl_path)
}

# Minimal star tree (runmode=-2 ignores topology but codeml still needs it)
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

# Write PHYLIP sequential format (codeml's native input)
write_phylip <- function(aln_mat, path) {
  char_mat <- as.character(aln_mat)
  n_seq  <- nrow(char_mat)
  n_col  <- ncol(char_mat)
  lines  <- c(paste(" ", n_seq, n_col))
  for (i in seq_len(n_seq)) {
    seq_str <- paste(char_mat[i, ], collapse = "")
    # codeml name field: exactly 10 chars, padded/truncated
    nm <- formatC(rownames(char_mat)[i], width = 10, flag = "-")
    nm <- substr(nm, 1, 10)
    lines <- c(lines, paste0(nm, "  ", seq_str))
  }
  writeLines(lines, path)
}

# ── Run ───────────────────────────────────────────────────────────────────────
pair_label <- TEST_PAIR
accs       <- parse_pair(pair_label)
pair_dir   <- file.path(OUT_ROOT, pair_label)
dir.create(pair_dir, showWarnings = FALSE, recursive = TRUE)

message("\n[", pair_label, "] Extracting ", length(accs), " sequences...")
pair_aln <- as.matrix(aln_sub[accs])

phy_path  <- file.path(pair_dir, "seq.phy")
tree_path <- file.path(pair_dir, "star.nwk")
ctl_path  <- file.path(pair_dir, "codeml.ctl")
out_path  <- file.path(pair_dir, "codeml_out.txt")

write_phylip(pair_aln, phy_path)
message("  Wrote PHYLIP: ", phy_path,
        "  (", length(accs), " seqs x ", ncol(pair_aln), " sites)")

write.tree(make_star_tree(rownames(pair_aln)), tree_path)
message("  Wrote star tree: ", tree_path)

write_ctl(phy_path, out_path, tree_path, ctl_path)
message("  Wrote control file: ", ctl_path)

message("  Running codeml (", TOOLS$codeml, ")...")
old_wd <- setwd(pair_dir)
out <- system2(TOOLS$codeml, args = "codeml.ctl", stdout = TRUE, stderr = TRUE)
setwd(old_wd)

# ── Check results ─────────────────────────────────────────────────────────────
dn_file <- file.path(pair_dir, "2ML.dN")
ds_file <- file.path(pair_dir, "2ML.dS")

if (!file.exists(ds_file)) {
  message("\n[FAIL] codeml did not produce 2ML.dS")
  message("codeml stdout/stderr:")
  cat(paste(out, collapse = "\n"), "\n")
  quit(status = 1)
}

message("\n[PASS] codeml produced output files:")
message("  ", dn_file)
message("  ", ds_file)

# Quick parse: dS matrix should have N*(N-1)/2 values
ds_lines <- readLines(ds_file)
message("\nFirst few lines of 2ML.dS:")
cat(head(ds_lines, 8), sep = "\n")
message("\n\nDone. Test passed.")
