# scripts/02_phylogenetics/01_maple_tree.R
#
# PURPOSE: Infer an initial phylogenetic tree with MAPLE, then prune outlier
#          taxa in two rounds based on branch length thresholds.
#
# MAPLE requires a custom input format. This script:
#   1. Converts the aligned FASTA → MAPLE format via createMapleFile.py
#   2. Runs MAPLEv0.6.8.py to infer the tree
#   3. Prunes outlier tips (very short / very long branches)
#
# INPUT:
#   data/alignments/{gene}_aligned.fasta
#
# OUTPUT:
#   data/trees/raw/maple/{gene}_tree.nwk         — raw MAPLE tree (Newick)
#   data/trees/pruned/{gene}_pruned.nwk          — pruned tree (Newick)
#   data/trees/pruned/{gene}_accessions.txt      — accessions retained
#   figures/trees/branch_hist.pdf                — raw branch length histograms
#   figures/trees/pruned_trees.pdf               — pruned tree plots

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
  library(ggplot2)
})

source(here("config.R"))
source(here("R", "tree_utils.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else ALL_GENES

maple_dir      <- TOOLS$maple
create_maple   <- file.path(maple_dir, "createMapleFile.py")
maple_script   <- file.path(maple_dir, "MAPLEv0.6.8.py")

raw_trees <- list()   # tree objects keyed by gene, for histogram + pruning

for (gene in genes) {
  message("\n[", gene, "] MAPLE tree inference...")

  fasta_in  <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  tree_out  <- file.path(PATHS$trees_maple, paste0(gene, "_tree.nwk"))
  acc_out   <- file.path(PATHS$trees_maple, paste0(gene, "_accessions_final.txt"))

  dir.create(PATHS$trees_maple, recursive = TRUE, showWarnings = FALSE)

  if (file.exists(tree_out)) {
    message("  Tree already exists, skipping MAPLE: ", tree_out)
    raw_trees[[gene]] <- read.tree(tree_out)
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
        fasta_in <- file.path(PATHS$trees_maple, paste0(gene, "_pdm09.fasta"))
        write.FASTA(seqs, fasta_in)
      }
    }
  }

  # ── Step 1: convert FASTA → MAPLE format ─────────────────────────────────
  # createMapleFile.py expects --path as the working dir,
  # with --fasta and --output as filenames relative to that path.
  fasta_copy <- file.path(PATHS$trees_maple, paste0(gene, "_for_maple.fasta"))
  file.copy(fasta_in, fasta_copy, overwrite = TRUE)

  maple_input_name <- paste0(gene, "_maple_input.txt")
  message("  Converting to MAPLE format...")

  ret <- system2("python3", args = c(
    create_maple,
    "--path", paste0(PATHS$trees_maple, "/"),
    "--fasta", paste0(gene, "_for_maple.fasta"),
    "--output", maple_input_name,
    "--overwrite"
  ))
  if (ret != 0) stop("createMapleFile.py failed for ", gene)
  unlink(fasta_copy)

  # ── Step 2: run MAPLE ────────────────────────────────────────────────────
  maple_input  <- file.path(PATHS$trees_maple, maple_input_name)
  maple_output <- file.path(PATHS$trees_maple, paste0(gene, "_maple"))
  message("  Running MAPLE...")

  ret <- system2("python3", args = c(
    maple_script,
    "--input", maple_input,
    "--output", maple_output,
    "--overwrite"
  ))
  if (ret != 0) stop("MAPLE failed for ", gene)

  # MAPLE outputs {output}_tree.tree
  maple_tree_file <- paste0(maple_output, "_tree.tree")
  if (!file.exists(maple_tree_file)) {
    stop("MAPLE output tree not found: ", maple_tree_file)
  }
  message("  MAPLE tree: ", basename(maple_tree_file))

  tree <- read.tree(maple_tree_file)
  message("  Raw tree: ", Ntip(tree), " tips")

  write.tree(tree, tree_out)
  writeLines(tree$tip.label, acc_out)
  message("  Tree:        ", tree_out)
  message("  Accessions:  ", acc_out)

  # ── Plot raw tree topology ──────────────────────────────────────────────
  fig_dir <- here("figures", "trees")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pdf_path <- file.path(fig_dir, paste0(gene, "_raw.pdf"))

  pdf(pdf_path, width = 10, height = 10)
  plot(tree, show.tip.label = FALSE, no.margin = TRUE,
       edge.width = 0.3, type = "phylogram",
       main = paste0(gene, " — raw MAPLE tree (", Ntip(tree), " tips)"))
  add.scale.bar(cex = 0.7)
  dev.off()
  message("  Plot:        ", pdf_path)

  raw_trees[[gene]] <- tree

  # Clean up intermediate files
  unlink(maple_input)
  maple_temps <- list.files(
    PATHS$trees_maple,
    pattern = paste0("^", gene, "_maple"),
    full.names = TRUE
  )
  unlink(maple_temps)
}

# ── Branch length histograms (raw trees) ─────────────────────────────────────
if (length(raw_trees) > 0) {
  fig_dir  <- here("figures", "trees")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  hist_pdf <- file.path(fig_dir, "branch_hist.pdf")

  pdf(hist_pdf, width = 8, height = 5)
  for (gene in names(raw_trees)) {
    bl <- raw_trees[[gene]]$edge.length
    p <- ggplot(data.frame(bl = bl), aes(x = bl)) +
      geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
      labs(title = paste0(gene, " — raw branch lengths"),
           x = "Branch length", y = "Count") +
      theme_bw()
    print(p)
  }
  dev.off()
  message("Branch length histograms: ", hist_pdf)
}

# ── Pruning ───────────────────────────────────────────────────────────────────
dir.create(PATHS$trees_pruned, recursive = TRUE, showWarnings = FALSE)
fig_dir <- here("figures", "trees")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

pruned_pdf <- file.path(fig_dir, "pruned_trees.pdf")
pdf(pruned_pdf, width = 10, height = 10)

for (gene in names(raw_trees)) {
  message("\n[", gene, "] Pruning...")
  tree <- raw_trees[[gene]]

  # Round 1: remove very short branches (likely duplicates/identical sequences)
  tree <- prune_outliers(tree,
                         min_bl = PARAMS$prune_min_bl,
                         max_bl = Inf,
                         round  = 1)

  # Round 2: remove very long branches (likely sequencing errors/mis-assignments)
  tree <- prune_outliers(tree,
                         min_bl = 0,
                         max_bl = PARAMS$prune_max_bl,
                         round  = 2)

  message("  Pruned tree: ", Ntip(tree), " tips")

  pruned_out <- file.path(PATHS$trees_pruned, paste0(gene, "_pruned.nwk"))
  acc_out    <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
  write.tree(tree, pruned_out)
  writeLines(tree$tip.label, acc_out)
  message("  Pruned tree: ", pruned_out)

  plot(tree, show.tip.label = FALSE, no.margin = TRUE,
       edge.width = 0.3, type = "phylogram",
       main = paste0(gene, " — pruned MAPLE tree (", Ntip(tree), " tips)"))
  add.scale.bar(cex = 0.7)
}

dev.off()
message("Pruned tree plots: ", pruned_pdf)

message("\nDone.")
