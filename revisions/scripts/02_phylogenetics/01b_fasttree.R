# scripts/02_phylogenetics/01b_fasttree.R
#
# PURPOSE: Infer a phylogenetic tree with FastTree (GTR+CAT) as an
#          alternative to MAPLE for divergent datasets like influenza HA/NA.
#
# FastTree is installed in the climflu conda env. Activate it before running:
#   conda activate climflu
#   Rscript scripts/02_phylogenetics/01b_fasttree.R
#
# INPUT:
#   data/alignments/{gene}_aligned.fasta
#
# OUTPUT:
#   data/trees/raw/fasttree/{gene}_fasttree.nwk             — tree (Newick)
#   data/trees/raw/fasttree/{gene}_fasttree_accessions.txt  — tip labels
#   figures/trees/{gene}_fasttree_raw.pdf                   — topology plot

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "tree_utils.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else ALL_GENES

# FastTree binary — use FastTreeMP (OpenMP) if available for threading
fasttree_bin <- Sys.which("FastTreeMP")
if (nchar(fasttree_bin) == 0) fasttree_bin <- Sys.which("FastTree")
if (nchar(fasttree_bin) == 0) fasttree_bin <- Sys.which("fasttree")
if (nchar(fasttree_bin) == 0) {
  stop("FastTree not found. Install via: mamba install -c bioconda fasttree")
}
message("Using: ", fasttree_bin)

for (gene in genes) {
  message("\n[", gene, "] FastTree inference...")

  fasta_in <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  dir.create(PATHS$trees_fast, recursive = TRUE, showWarnings = FALSE)
  tree_out <- file.path(PATHS$trees_fast, paste0(gene, "_fasttree.nwk"))
  acc_out  <- file.path(PATHS$trees_fast, paste0(gene, "_fasttree_accessions.txt"))

  if (file.exists(tree_out)) {
    message("  Tree already exists, skipping: ", tree_out)
    next
  }

  if (!file.exists(fasta_in)) {
    warning("Alignment not found, skipping: ", fasta_in)
    next
  }

  # ── Filter: remove seasonal H1N1 strains (keep pdm09 only) ───────────────
  h1n1_genes <- c("H1", "N1", "PA_H1N1")
  if (gene %in% h1n1_genes) {
    climate_csv <- file.path(PATHS$climate, paste0(gene, "_climate.csv"))
    if (file.exists(climate_csv)) {
      climate <- read_csv(climate_csv, show_col_types = FALSE)
      pdm_acc <- climate$accession[climate$lineage == "pdm09" &
                                     !is.na(climate$lineage)]
      seqs <- read.FASTA(fasta_in)
      keep <- intersect(names(seqs), pdm_acc)
      n_dropped <- length(seqs) - length(keep)
      if (n_dropped > 0) {
        message("  Removed ", n_dropped, " seasonal H1N1 strains ",
                "(keeping ", length(keep), " pdm09)")
        seqs <- seqs[keep]
        fasta_in <- file.path(PATHS$trees_fast, paste0(gene, "_pdm09.fasta"))
        write.FASTA(seqs, fasta_in)
      }
    }
  }

  # ── Run FastTree ─────────────────────────────────────────────────────────
  # -gtr: generalised time-reversible model
  # -nt:  nucleotide alignment
  # -fastest: speed up for large alignments (NNI + SPR)
  log_file <- file.path(PATHS$trees_fast, paste0(gene, "_fasttree.log"))
  message("  Running FastTree (GTR+CAT)...")
  ret <- system2(fasttree_bin,
                 args = c("-gtr", "-nt", "-fastest", fasta_in),
                 stdout = tree_out,
                 stderr = log_file)

  if (ret != 0 || !file.exists(tree_out) || file.size(tree_out) == 0) {
    warning("FastTree failed for ", gene, " (see ", log_file, ")")
    next
  }

  tree <- read.tree(tree_out)
  message("  Tree: ", Ntip(tree), " tips")

  writeLines(tree$tip.label, acc_out)

  # ── Plot raw tree topology ──────────────────────────────────────────────
  fig_dir <- here("revisions", "figures", "trees")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(fig_dir, paste0(gene, "_fasttree_raw.pdf"))

  pdf(pdf_path, width = 10, height = max(8, Ntip(tree) * 0.01))
  plot(tree, show.tip.label = FALSE, no.margin = TRUE,
       edge.width = 0.3, type = "phylogram",
       main = paste0(gene, " — raw FastTree (",
                     Ntip(tree), " tips)"))
  add.scale.bar(cex = 0.7)
  dev.off()
  message("  Plot: ", pdf_path)
}

message("\nDone.")
