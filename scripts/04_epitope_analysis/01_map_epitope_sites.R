# scripts/04_epitope_analysis/01_map_epitope_sites.R
#
# PURPOSE: Map known antigenic site positions onto the MAFFT alignment columns
#          for H1 and H3 HA using canonical reference sequences.
#          The references are added to a TEMP copy of the alignment via
#          mafft --add --keeplength — the original alignment is NEVER modified.
#
# REFERENCES (fetched via ape::read.GenBank):
#   H3: A/Aichi/2/1968    CY121680  (canonical H3 numbering)
#   H1: A/Puerto Rico/8/1934  CY009444  (canonical H1 numbering)
#
# ANTIGENIC SITES:
#   H3: sites A–E  (Wilson & Cox 1990, mature H3 numbering)
#   H1: Sa, Sb, Ca1, Ca2, Cb  (Caton et al. 1982, mature H1 numbering)
#
# ANTIGENIC RIDGE (subset of epitope most responsible for antigenic transitions):
#   H3: 7 key positions from Koel et al. 2013 (Science 342:976-979)
#   H1: 7 key Sa/Sb ridge positions from Xu et al. 2010 + Caton et al. 1982
#
# OUTPUT:
#   data/epitope/{gene}_epitope_cols.txt      alignment column indices (1-based)
#   data/epitope/{gene}_nonepitope_cols.txt
#   data/epitope/{gene}_ridge_cols.txt        antigenic ridge only
#   data/epitope/{gene}_nonridge_cols.txt     all non-ridge columns
#   data/epitope/{gene}_site_map.csv          per-column annotation table

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
})

source(here("config.R"))

out_dir <- here("data", "epitope")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Antigenic site definitions (mature protein numbering) ─────────────────────
# H3: Wilson & Cox 1990 / Wiley & Skehel 1987
antigenic_sites_h3 <- list(
  A = c(122, 124, 126, 130, 131, 132, 133, 135, 137, 138,
        140, 142, 143, 144, 145, 146, 150, 152, 168),
  B = c(128, 129, 155, 156, 157, 158, 159, 160, 163, 165,
        186, 187, 188, 189, 190, 192, 193, 194, 196, 197, 198),
  C = c(44, 45, 46, 47, 48, 50, 51, 53, 54,
        273, 275, 276, 278, 279, 280, 294, 297, 299),
  D = c(96, 102, 103, 117, 121, 167, 170, 171, 172, 173,
        174, 175, 176, 177, 179, 182, 185, 222, 224, 225,
        228, 233, 234),
  E = c(57, 59, 62, 63, 67, 75, 78, 80, 81, 82, 83,
        86, 87, 88, 91, 92, 94, 109, 260, 261, 262, 265)
)

# H1: Caton et al. 1982 + Gerhard et al. 1981 (mature H1 numbering)
antigenic_sites_h1 <- list(
  Sa  = c(124, 125, 126, 127, 128, 129, 130, 131, 132, 133,
          134, 135, 136, 137, 138, 139, 140, 141),
  Sb  = c(156, 157, 158, 159, 160, 161, 162, 163, 164,
          165, 166, 167, 168),
  Ca1 = c(166, 167, 168, 169, 170, 171, 172, 173,
          174, 175, 176, 177),
  Ca2 = c(203, 204, 205, 206, 207, 208, 209, 210, 211),
  Cb  = c(73, 74, 75, 76, 77, 78, 79, 80, 81)
)

# Antigenic ridge: minimal set of positions driving most antigenic transitions.
# H3: Koel et al. 2013 (Science 342:976-979) — 7 positions.
# H1: key Sa/Sb ridge residues from Xu et al. 2010 + Caton et al. 1982.
antigenic_ridge_h3 <- c(145, 155, 156, 158, 159, 189, 193)
antigenic_ridge_h1 <- c(133, 137, 138, 140, 141, 156, 157)

# Signal peptide length (aa) and HA1/HA2 cleavage site (last HA1 residue,
# mature protein numbering)
gene_meta <- list(
  H3 = list(accession   = "CY121680",
             ref_name    = "ref_H3_Aichi68",
             signal_aa   = 16L,
             ha1_last_aa = 329L,
             sites       = antigenic_sites_h3,
             ridge       = antigenic_ridge_h3),
  H1 = list(accession   = "CY009444",
             ref_name    = "ref_H1_PR8",
             signal_aa   = 17L,
             ha1_last_aa = 327L,
             sites       = antigenic_sites_h1,
             ridge       = antigenic_ridge_h1)
)

# ── Helper: fetch reference from GenBank ──────────────────────────────────────
fetch_reference <- function(accession, ref_name, cache_dir) {
  cache_path <- file.path(cache_dir, paste0(ref_name, ".fasta"))
  if (file.exists(cache_path)) {
    message("  Using cached reference: ", cache_path)
    return(cache_path)
  }
  message("  Fetching ", accession, " from NCBI GenBank...")
  seq <- read.GenBank(accession, as.character = TRUE)
  seq_str <- paste(seq[[1]], collapse = "")
  writeLines(c(paste0(">", ref_name), seq_str), cache_path)
  message("  Written: ", cache_path)
  cache_path
}

# ── Helper: add reference to alignment via mafft --add (temp file only) ───────
add_reference <- function(aln_path, ref_fasta, mafft_bin, threads = 1L) {
  tmp_out <- tempfile(fileext = ".fasta")
  ret <- system2(mafft_bin,
                 args   = c("--quiet", "--keeplength", "--thread", threads,
                            "--add", ref_fasta, aln_path),
                 stdout = tmp_out, stderr = FALSE)
  if (ret != 0 || !file.exists(tmp_out))
    stop("mafft --add failed (return code ", ret, ")")
  tmp_out
}

# ── Helper: build codon→column lookup from reference row ──────────────────────
# Returns a named integer matrix: rownames = mature aa position (may be
# negative for signal peptide), 3 columns = alignment column indices for
# each nt in that codon.
build_codon_map <- function(aln, ref_name, signal_aa) {
  if (!ref_name %in% rownames(aln))
    stop("Reference '", ref_name, "' not found in extended alignment")
  ref_chars  <- aln[ref_name, ]
  non_gap    <- which(!(ref_chars %in% c("-", "n", "N", "?")))
  n_codons   <- floor(length(non_gap) / 3)
  codon_mat  <- matrix(non_gap[seq_len(n_codons * 3)],
                       ncol = 3, byrow = TRUE)
  # Row = precursor aa position (1-based); mature = row - signal_aa
  mature_pos <- seq_len(nrow(codon_mat)) - signal_aa
  rownames(codon_mat) <- as.character(mature_pos)
  codon_mat
}

# ── Helper: resolve site positions to alignment columns ───────────────────────
sites_to_cols <- function(sites_list, codon_map) {
  cols <- lapply(names(sites_list), function(site_name) {
    positions <- sites_list[[site_name]]
    cols_for_site <- unlist(lapply(as.character(positions), function(p) {
      if (p %in% rownames(codon_map)) codon_map[p, ] else integer(0)
    }))
    data.frame(site = site_name,
               residue_pos = rep(as.integer(positions),
                                 each = 3)[seq_along(cols_for_site)],
               alignment_col = cols_for_site)
  })
  bind_rows(cols)
}

# ── Main loop ─────────────────────────────────────────────────────────────────
for (gene in c("H1", "H3")) {
  message("\n[", gene, "] Mapping epitope sites...")

  meta     <- gene_meta[[gene]]
  aln_path <- file.path(PATHS$alignments, paste0(gene, "_aligned.fasta"))
  if (!file.exists(aln_path))
    stop("Alignment not found: ", aln_path)

  # Fetch/cache reference
  ref_fasta  <- fetch_reference(meta$accession, meta$ref_name, out_dir)

  # Add reference to alignment (temp copy — original untouched)
  message("  Running mafft --add --keeplength (temp file)...")
  tmp_ext <- add_reference(aln_path, ref_fasta,
                           TOOLS$mafft, threads = PARAMS$mafft_threads)

  # Read extended alignment
  ext_aln <- read.dna(tmp_ext, format = "fasta", as.character = TRUE,
                      as.matrix = TRUE)
  n_orig_cols <- ncol(read.dna(aln_path, format = "fasta",
                               as.character = TRUE, as.matrix = TRUE))
  n_ext_cols  <- ncol(ext_aln)
  message("  Original alignment columns: ", n_orig_cols,
          "  Extended: ", n_ext_cols)
  if (n_ext_cols != n_orig_cols)
    warning("  Column count changed after mafft --add — ",
            "reference may have insertions not in alignment. ",
            "Using only original-alignment columns.")

  # Build codon→column map from reference row
  codon_map <- build_codon_map(ext_aln, meta$ref_name, meta$signal_aa)
  message("  Reference codons mapped: ", nrow(codon_map),
          "  (signal peptide: ", meta$signal_aa, " aa)")

  # Resolve epitope site positions to alignment columns
  epitope_df <- sites_to_cols(meta$sites, codon_map) |>
    filter(alignment_col <= n_orig_cols)   # discard any ref-only columns
  message("  Epitope alignment columns: ", n_distinct(epitope_df$alignment_col))

  epitope_cols    <- sort(unique(epitope_df$alignment_col))
  nonepitope_cols <- setdiff(seq_len(n_orig_cols), epitope_cols)
  message("  Non-epitope columns: ", length(nonepitope_cols))

  # Antigenic ridge columns
  ridge_df <- data.frame(
    residue_pos   = rep(meta$ridge, each = 3),
    alignment_col = unlist(lapply(as.character(meta$ridge), function(p) {
      if (p %in% rownames(codon_map)) codon_map[p, ] else integer(0)
    }))
  ) |> filter(alignment_col <= n_orig_cols)
  ridge_cols    <- sort(unique(ridge_df$alignment_col))
  nonridge_cols <- setdiff(seq_len(n_orig_cols), ridge_cols)
  message("  Ridge columns: ", length(ridge_cols),
          "  (", length(meta$ridge), " positions)")

  # ── Build site map table ───────────────────────────────────────────────────
  all_mature <- as.integer(rownames(codon_map))
  site_map <- lapply(seq_len(nrow(codon_map)), function(i) {
    mature  <- all_mature[i]
    cols    <- codon_map[i, ]
    cols    <- cols[cols <= n_orig_cols]  # drop any ref-only extras
    if (length(cols) == 0) return(NULL)
    # Region
    region <- if (mature <= 0) {
      "signal"
    } else if (mature <= meta$ha1_last_aa) {
      "HA1"
    } else {
      "HA2"
    }
    # Site name (if epitope); flag ridge positions
    site_hit <- names(Filter(function(s) mature %in% s, meta$sites))
    site_name <- if (length(site_hit) > 0) site_hit[1] else NA_character_
    if (!is.na(site_name)) region <- paste0("epitope_", site_name)
    if (mature %in% meta$ridge) region <- paste0(region, "_ridge")

    data.frame(
      alignment_col = cols,
      nt_start      = cols,
      nt_end        = cols,
      residue_pos   = mature,
      region        = region,
      site          = site_name,
      stringsAsFactors = FALSE
    )
  }) |> bind_rows()

  # ── Write outputs ──────────────────────────────────────────────────────────
  writeLines(as.character(epitope_cols),
             file.path(out_dir, paste0(gene, "_epitope_cols.txt")))
  writeLines(as.character(nonepitope_cols),
             file.path(out_dir, paste0(gene, "_nonepitope_cols.txt")))
  writeLines(as.character(ridge_cols),
             file.path(out_dir, paste0(gene, "_ridge_cols.txt")))
  writeLines(as.character(nonridge_cols),
             file.path(out_dir, paste0(gene, "_nonridge_cols.txt")))
  write_csv(site_map,
            file.path(out_dir, paste0(gene, "_site_map.csv")))
  message("  Written: ", gene,
          "_epitope_cols.txt / _nonepitope_cols.txt / ",
          "_ridge_cols.txt / _nonridge_cols.txt / _site_map.csv")

  unlink(tmp_ext)   # clean up temp extended alignment
}

message("\nDone.")
