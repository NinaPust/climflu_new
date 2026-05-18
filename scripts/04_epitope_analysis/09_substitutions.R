# scripts/04_epitope_analysis/09_substitutions.R
#
# PURPOSE: Count NS and S substitutions on terminal branches via ancestral
#          parsimony (phangorn::ancestral.pars) for any combination of:
#            ridge   — 7 antigenic-ridge codons
#            epitope — full epitope alignment (~103 codons)
#            full    — whole HA1+HA2 gene
#            HA2     — HA2 subunit (conserved; negative control)
#
#          Results are appended to data/rates/{gene}_rates_epitope.csv.
#          Mutation records (individual NS events with residue positions)
#          are saved for ridge and epitope regions.
#
# INPUT:
#   data/epitope/{gene}_{region}_aligned.fasta
#   data/alignments/{gene}_aligned.fasta          (full only)
#   data/trees/pruned/{gene}_accessions.txt        (full only)
#   data/trees/timed/{gene}_timetree.nwk
#   data/rates/{gene}_rates_epitope.csv
#   data/epitope/{gene}_{region}_cols.txt          (epitope/full)
#   data/epitope/{gene}_site_map.csv               (epitope/full)
#
# OUTPUT:
#   data/rates/{gene}_rates_epitope.csv            (new columns appended)
#   data/epitope/{gene}_{region}_mutations.csv     (ridge and epitope)
#
# FLAGS:
#   --gene X      run a single gene (default: H1 H3)
#   --region X    ridge | epitope | full  (default: ridge epitope full)
#   --cores N     parallel cores for epitope/full (default: 8)
#   --sequential  force 1 core

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(phangorn)
  library(Biostrings)
  library(parallel)
  library(dplyr)
  library(readr)
  library(tibble)
})

source(here("config.R"))

args    <- commandArgs(trailingOnly = TRUE)
genes   <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
regions <- if ("--region" %in% args) args[which(args == "--region") + 1] else c("ridge", "epitope", "full", "HA2")
n_cores <- if ("--cores"  %in% args) as.integer(args[which(args == "--cores") + 1]) else
           if ("--sequential" %in% args) 1L else 8L

# ── Codon helpers ─────────────────────────────────────────────────────────────

genetic_code <- as.list(Biostrings::GENETIC_CODE)

codon_to_aa <- function(codon) {
  codon <- toupper(codon)
  if (any(c("N", "-", "?") %in% strsplit(codon, "")[[1]])) return(NA_character_)
  aa <- genetic_code[[codon]]
  if (is.null(aa)) NA_character_ else aa
}

classify_change <- function(c1, c2) {
  if (length(c1) == 0 || length(c2) == 0) return(NULL)
  if (is.na(c1) || is.na(c2))             return(NULL)
  if (c1 == c2) return(list(S = 0L, NS = 0L))
  aa1 <- codon_to_aa(c1); aa2 <- codon_to_aa(c2)
  if (is.na(aa1) || is.na(aa2))           return(NULL)
  if (aa1 == "*" || aa2 == "*")           return(NULL)
  if (aa1 == aa2) list(S = 1L, NS = 0L) else list(S = 0L, NS = 1L)
}

# Process one codon position; returns list(NS, S, mutations data.frame).
# pos: mature residue position (integer, may be NA if unknown).
process_codon <- function(ci, aln, tree, tip_parent, n_tips, pos = NA_integer_) {
  cols       <- ((ci - 1L) * 3L + 1L):(ci * 3L)
  codon_seqs <- apply(aln[, cols, drop = FALSE], 1,
                      function(x) paste(toupper(x), collapse = ""))
  all_codons   <- sort(unique(codon_seqs))
  phyd         <- phyDat(as.list(codon_seqs), type = "USER", levels = all_codons)
  codon_levels <- attr(phyd, "levels")
  anc          <- ancestral.pars(tree, phyd, type = "MPR")

  parent_mat    <- do.call(rbind, lapply(tip_parent, function(i) anc[[i]]))
  parent_codons <- codon_levels[apply(parent_mat, 1, function(x) which(x > 0)[1])]
  tip_codons    <- codon_seqs[tree$tip.label]

  ns       <- integer(n_tips)
  s        <- integer(n_tips)
  mut_rows <- vector("list", n_tips)

  for (ti in seq_len(n_tips)) {
    ch <- classify_change(parent_codons[ti], tip_codons[ti])
    if (!is.null(ch)) {
      ns[ti] <- ch$NS
      s[ti]  <- ch$S
      if (ch$NS == 1L) {
        pa <- codon_to_aa(parent_codons[ti])
        ta <- codon_to_aa(tip_codons[ti])
        mut_rows[[ti]] <- data.frame(
          accession    = tree$tip.label[ti],
          codon_idx    = ci,
          residue_pos  = pos,
          parent_codon = parent_codons[ti],
          tip_codon    = tip_codons[ti],
          parent_aa    = pa,
          tip_aa       = ta,
          stringsAsFactors = FALSE
        )
      }
    }
  }
  list(NS = ns, S = s, mutations = do.call(rbind, mut_rows))
}

# Ridge positions per gene (must match 01_map_epitope_sites.R)
ridge_positions <- list(
  H3 = c(145, 155, 156, 158, 159, 189, 193),
  H1 = c(133, 137, 138, 140, 141, 156, 157)
)

# ── Main loop ─────────────────────────────────────────────────────────────────

for (gene in genes) {
  for (region in regions) {
    tag <- paste0("[", gene, " ", region, "]")
    message("\n", tag, " Counting substitutions via ancestral parsimony...")

    # ── Paths ─────────────────────────────────────────────────────────────────
    aln_path <- if (region == "full") {
      file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
    } else {
      file.path(here("data", "epitope"),
                paste0(gene, "_", region, "_aligned.fasta"))
    }
    tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
    rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))

    for (p in c(aln_path, tree_path, rates_path)) {
      if (!file.exists(p)) stop("Required file not found: ", p)
    }

    # ── Load data ─────────────────────────────────────────────────────────────
    tree <- read.tree(tree_path)
    aln  <- read.dna(aln_path, format = "fasta",
                     as.character = TRUE, as.matrix = TRUE)

    if (region == "full") {
      acc_path <- file.path(PATHS$trees_pruned, paste0(gene, "_accessions.txt"))
      if (file.exists(acc_path)) {
        pruned <- readLines(acc_path)
        pruned <- pruned[nzchar(trimws(pruned))]
        aln    <- aln[intersect(pruned, rownames(aln)), , drop = FALSE]
      }
    }

    common <- intersect(tree$tip.label, rownames(aln))
    tree   <- keep.tip(tree, common)
    aln    <- aln[common, , drop = FALSE]

    n_tips       <- Ntip(tree)
    n_codons     <- ncol(aln) / 3L
    tip_edge_idx <- match(seq_len(n_tips), tree$edge[, 2])
    tip_bl       <- setNames(tree$edge.length[tip_edge_idx], tree$tip.label)
    tip_parent   <- tree$edge[tip_edge_idx, 1]

    message("  Tips: ", n_tips, "  Codons: ", n_codons)

    # ── Build codon → residue position map ────────────────────────────────────
    # Ridge: positions from hardcoded list.
    # Epitope: derive from cols file + site_map.
    # Full: positions not needed (mutations not saved for full gene).
    codon_pos <- if (region == "ridge") {
      ridge_positions[[gene]]
    } else if (region == "full") {
      rep(NA_integer_, n_codons)
    } else {
      cols_file     <- file.path(here("data", "epitope"),
                                 paste0(gene, "_", region, "_cols.txt"))
      site_map_path <- file.path(here("data", "epitope"),
                                 paste0(gene, "_site_map.csv"))
      if (!file.exists(cols_file) || !file.exists(site_map_path)) {
        message("  site_map/cols not found, residue_pos will be NA")
        rep(NA_integer_, n_codons)
      } else {
        aln_cols     <- as.integer(readLines(cols_file))
        codon_starts <- aln_cols[seq(1, length(aln_cols), by = 3)]
        sm <- read_csv(site_map_path, show_col_types = FALSE) |>
          select(alignment_col, residue_pos)
        tibble(alignment_col = codon_starts) |>
          left_join(sm, by = "alignment_col") |>
          pull(residue_pos)
      }
    }

    # ── Run ancestral reconstruction ──────────────────────────────────────────
    cores_to_use <- if (region == "ridge") 1L else n_cores
    message("  Running ", n_codons, " codons on ", cores_to_use, " core(s)...")
    t0 <- proc.time()

    codon_results <- mclapply(
      seq_len(n_codons),
      function(ci) process_codon(ci, aln, tree, tip_parent, n_tips,
                                 codon_pos[ci]),
      mc.cores = cores_to_use
    )

    elapsed <- round((proc.time() - t0)[["elapsed"]])
    message("  Done in ", elapsed, "s")

    # ── Aggregate counts ──────────────────────────────────────────────────────
    ns_count <- Reduce("+", lapply(codon_results, `[[`, "NS"))
    s_count  <- Reduce("+", lapply(codon_results, `[[`, "S"))
    names(ns_count) <- tree$tip.label
    names(s_count)  <- tree$tip.label

    valid_bl <- tip_bl > 0
    ns_rate  <- ifelse(valid_bl, ns_count / tip_bl, NA_real_)
    s_rate   <- ifelse(valid_bl, s_count  / tip_bl, NA_real_)

    message("  Tips with NS > 0: ", sum(ns_count > 0),
            "  with S > 0: ", sum(s_count > 0))

    # Column prefix: "ridge", "epitope", "full"
    pfx <- region
    new_cols <- tibble(
      accession                           = names(ns_count),
      !!paste0(pfx, "_ns_count")       := ns_count,
      !!paste0(pfx, "_s_count")        := s_count,
      !!paste0(pfx, "_ns_rate")        := ns_rate,
      !!paste0(pfx, "_s_rate")         := s_rate,
      !!paste0("log_", pfx, "_ns_rate") :=
        ifelse(ns_count > 0, log(ns_rate), NA_real_),
      !!paste0("log_", pfx, "_s_rate") :=
        ifelse(s_count > 0, log(s_rate), NA_real_)
    )

    base   <- read_csv(rates_path, show_col_types = FALSE)
    base   <- base |> select(-any_of(names(new_cols)[-1]))
    df_out <- left_join(base, new_cols, by = "accession")
    write_csv(df_out, rates_path)
    message("  Updated: ", rates_path)

    # ── Save mutations CSV (ridge and epitope only) ───────────────────────────
    if (region %in% c("ridge", "epitope")) {
      mut_df <- do.call(rbind, lapply(codon_results, `[[`, "mutations"))
      if (!is.null(mut_df) && nrow(mut_df) > 0) {
        mut_df <- as_tibble(mut_df) |>
          mutate(
            gene   = gene,
            region = region,
            change = paste0(parent_aa, residue_pos, tip_aa)
          ) |>
          select(accession, gene, region, residue_pos,
                 parent_aa, tip_aa, parent_codon, tip_codon, change)
        mut_path <- file.path(here("data", "epitope"),
                              paste0(gene, "_", region, "_mutations.csv"))
        write_csv(mut_df, mut_path)
        message("  Mutations saved: ", nrow(mut_df), " NS events → ", mut_path)
      }
    }
  }
}

message("\nDone.")
