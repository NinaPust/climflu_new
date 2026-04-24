# scripts/02_phylogenetics/03_maple_r2.R
#
# PURPOSE: Build the canonical MAPLE tree by excluding both branch-length
#          outliers (from r1 pruning) and clock outliers (from 02_treetime_r1.R).
#          The resulting pruned r2 tree is the input for all downstream scripts.
#
#   1. Reads r1 pruned accessions and r1 TreeTime outliers.tsv
#   2. Subsets the alignment to (r1 accessions - clock outliers)
#   3. Runs MAPLE on the filtered alignment
#   4. Prunes with the same MAD-based method as 01_maple_r1.R
#
# INPUT:
#   data/trees/pruned/{gene}_accessions_r1.txt
#   data/trees/timed/{gene}_treetime_r1_out/outliers.tsv
#   data/trees/raw/maple/{gene}_pdm09_r1.fasta   (H1/N1/PA_H1N1 only)
#   data/alignments/{gene}_aligned.fasta         (H3/N2/PA_H3N2)
#
# OUTPUT:
#   data/trees/raw/maple/{gene}_tree_r2.nwk      — raw r2 MAPLE tree
#   data/trees/pruned/{gene}_pruned.nwk          — canonical pruned tree
#   data/trees/pruned/{gene}_accessions.txt      — accessions retained
#   figures/trees/branch_hist_r2.pdf             — r2 branch length histograms
#   figures/trees/pruned_trees_r2.pdf            — pruned r2 tree plots
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(readr)
  library(ggplot2)
})

source(here("config.R"))
source(here("R", "tree_utils.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args      <- commandArgs(trailingOnly = TRUE)
genes     <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}

maple_dir    <- TOOLS$maple
create_maple <- file.path(maple_dir, "createMapleFile.py")
maple_script <- file.path(maple_dir, "MAPLEv0.6.8.py")

r2_raw_trees <- list()

for (gene in genes) {
  message("\n[", gene, "] MAPLE r2 — clock-filtered alignment...")

  acc_r1_path    <- file.path(PATHS$trees_pruned,
                              paste0(gene, "_accessions_r1.txt"))
  outlier_path   <- file.path(PATHS$trees_timed,
                              paste0(gene, "_treetime_r1_out"), "outliers.tsv")
  tree_r2_out    <- file.path(PATHS$trees_maple, paste0(gene, "_tree_r2.nwk"))

  if (!file.exists(acc_r1_path)) {
    warning("R1 accessions not found, skipping: ", acc_r1_path)
    next
  }

  # ── Load r2 raw tree if it already exists ────────────────────────────────
  if (file.exists(tree_r2_out)) {
    message("  R2 tree already exists, skipping MAPLE: ", tree_r2_out)
    r2_raw_trees[[gene]] <- read.tree(tree_r2_out)
    next
  }

  # ── Build filtered accession list ────────────────────────────────────────
  r1_acc <- readLines(acc_r1_path)

  clock_outliers <- character(0)
  if (file.exists(outlier_path)) {
    ol <- read.delim(outlier_path, stringsAsFactors = FALSE)
    clock_outliers <- ol[[1]][!is.na(ol[[1]]) & ol[[1]] != ""]
    message("  Excluding ", length(clock_outliers), " clock outlier(s): ",
            paste(head(clock_outliers, 5), collapse = ", "),
            if (length(clock_outliers) > 5) "..." else "")
  } else {
    message("  No outliers.tsv found — using full r1 accession set")
  }

  keep_acc <- setdiff(r1_acc, clock_outliers)
  message("  Alignment subset: ", length(keep_acc), " sequences ",
          "(", length(clock_outliers), " clock outliers removed)")

  # ── Source alignment ──────────────────────────────────────────────────────
  pdm09_fasta <- file.path(PATHS$trees_maple,
                           paste0(gene, "_pdm09_r1.fasta"))
  src_fasta   <- if (file.exists(pdm09_fasta)) {
    pdm09_fasta
  } else {
    file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  }

  if (!file.exists(src_fasta)) {
    warning("Source alignment not found, skipping: ", src_fasta)
    next
  }

  seqs <- read.FASTA(src_fasta)
  seqs <- seqs[intersect(names(seqs), keep_acc)]

  r2_fasta <- file.path(PATHS$trees_maple,
                        paste0(gene, "_for_maple_r2.fasta"))
  write.FASTA(seqs, r2_fasta)

  # ── Convert to MAPLE format ───────────────────────────────────────────────
  maple_input_name <- paste0(gene, "_maple_input_r2.txt")
  message("  Converting to MAPLE format...")

  ret <- system2("python3", args = c(
    create_maple,
    "--path",   paste0(PATHS$trees_maple, "/"),
    "--fasta",  paste0(gene, "_for_maple_r2.fasta"),
    "--output", maple_input_name,
    "--overwrite"
  ))
  if (ret != 0) {
    warning("createMapleFile.py failed for ", gene)
    next
  }
  unlink(r2_fasta)

  # ── Run MAPLE ─────────────────────────────────────────────────────────────
  maple_input  <- file.path(PATHS$trees_maple, maple_input_name)
  maple_output <- file.path(PATHS$trees_maple, paste0(gene, "_maple_r2"))
  message("  Running MAPLE...")

  ret <- system2("python3", args = c(
    maple_script,
    "--input",  maple_input,
    "--output", maple_output,
    "--overwrite"
  ))
  if (ret != 0) {
    warning("MAPLE failed for ", gene)
    next
  }

  maple_tree_file <- paste0(maple_output, "_tree.tree")
  if (!file.exists(maple_tree_file)) {
    warning("MAPLE output not found: ", maple_tree_file)
    next
  }

  tree <- read.tree(maple_tree_file)
  message("  Raw r2 tree: ", Ntip(tree), " tips")

  write.tree(tree, tree_r2_out)
  message("  Written: ", tree_r2_out)

  # ── Plot raw r2 tree ──────────────────────────────────────────────────────
  fig_dir  <- here("figures", "trees")
  dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
  pdf(file.path(fig_dir, paste0(gene, "_raw_r2.pdf")), width = 10, height = 10)
  plot(tree, show.tip.label = FALSE, no.margin = TRUE,
       edge.width = 0.3, type = "phylogram",
       main = paste0(gene, " — raw MAPLE tree r2 (", Ntip(tree), " tips)"))
  add.scale.bar(cex = 0.7)
  dev.off()

  r2_raw_trees[[gene]] <- tree

  unlink(maple_input)
  unlink(list.files(PATHS$trees_maple,
                    pattern    = paste0("^", gene, "_maple_r2"),
                    full.names = TRUE))
}

# ── Branch length histograms (r2 raw trees) ───────────────────────────────────
if (length(r2_raw_trees) > 0) {
  fig_dir  <- here("figures", "trees")
  hist_pdf <- file.path(fig_dir, "branch_hist_r2.pdf")
  pdf(hist_pdf, width = 8, height = 5)
  for (gene in names(r2_raw_trees)) {
    bl <- r2_raw_trees[[gene]]$edge.length
    p <- ggplot(data.frame(bl = bl), aes(x = bl)) +
      geom_histogram(bins = 60, fill = "steelblue", colour = "white") +
      labs(title = paste0(gene, " — r2 raw branch lengths"),
           x = "Branch length", y = "Count") +
      theme_bw()
    print(p)
  }
  dev.off()
  message("R2 branch length histograms: ", hist_pdf)
}

# ── Prune r2 trees → canonical output ────────────────────────────────────────
dir.create(PATHS$trees_pruned, recursive = TRUE, showWarnings = FALSE)
fig_dir    <- here("figures", "trees")
pruned_pdf <- file.path(fig_dir, "pruned_trees_r2.pdf")
pdf(pruned_pdf, width = 10, height = 10)

for (gene in names(r2_raw_trees)) {
  message("\n[", gene, "] Pruning r2 tree...")
  tree <- r2_raw_trees[[gene]]

  tree <- prune_outliers(tree,
                         min_bl = PARAMS$prune_min_bl,
                         k_long = Inf,
                         round  = 1)
  tree <- prune_outliers(tree,
                         min_bl = 0,
                         k_long = PARAMS$prune_long_k,
                         round  = 2)

  message("  Pruned tree: ", Ntip(tree), " tips")

  pruned_out <- file.path(PATHS$trees_pruned, paste0(gene, "_pruned.nwk"))
  acc_out    <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
  write.tree(tree, pruned_out)
  writeLines(tree$tip.label, acc_out)
  message("  Written: ", pruned_out)

  plot(tree, show.tip.label = FALSE, no.margin = TRUE,
       edge.width = 0.3, type = "phylogram",
       main = paste0(gene, " — pruned r2 MAPLE tree (", Ntip(tree), " tips)"))
  add.scale.bar(cex = 0.7)
}

dev.off()
message("Pruned r2 tree plots: ", pruned_pdf)

message("\nDone.")
