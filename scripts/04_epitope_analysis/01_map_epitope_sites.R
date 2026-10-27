# scripts/04_epitope_analysis/01_map_epitope_sites.R
#
# PURPOSE: Map antigenic site positions onto MAFFT alignment columns for H1 and
#          H3 HA, using canonical reference sequences. Produces the column
#          definitions every downstream subregion analysis depends on.
#
# This is a rewrite. Three defects in the original version are fixed here, all
# of which silently corrupted the region definitions:
#
#   1. WRONG H3 REFERENCE. The H3 accession was CY121680, which returns an H1
#      HA, not A/Aichi/2/1968. Every H3 antigenic site was therefore mapped
#      through an H1 sequence. Symptoms: 156 of 1752 reference bases with no
#      alignment column, the fusion peptide unlocatable, and 15% of epitope
#      positions lost. Corrected accession: EF614248.1.
#      Guarded by verify_reference_subtype() below.
#
#   2. FRAME MIS-REGISTRATION. Codons were built from the reference's FIRST
#      base, assuming the FASTA starts at the CDS. Both references carry a 5'
#      leader (19 nt for H1), so every codon boundary was out of phase — for
#      every region, not just full/HA2. Translating H1 at the old mapping gave
#      40 internal stops; the corrected mapping gives the canonical signal
#      peptide (MKANLLVLLCALAAADA) and fusion peptide (GLFGAIAGFIEG).
#      The CDS is now located on the ORIGINAL reference FASTA, where both
#      references have clean full-length ORFs.
#
#   3. REFERENCE MUTILATION BY --keeplength. That flag forces the added
#      reference into the existing column count and silently DELETES bases
#      that would need a new column, costing H1 its start codon. Replaced with
#      --add --mapout, which makes MAFFT emit an authoritative base-to-column
#      map with unmappable bases explicitly marked.
#
# Reference caches are keyed by accession, so changing an accession forces a
# re-fetch rather than silently reusing a stale (or wrong) cached file.
#
# REFERENCES:
#   H3: A/Aichi/2/1968        EF614248.1  (canonical H3 numbering)
#   H1: A/Puerto Rico/8/1934  CY009444    (canonical H1 numbering)
#
# ANTIGENIC SITES:
#   H3: sites A-E  (Wilson & Cox 1990, mature H3 numbering)
#   H1: Sa, Sb, Ca1, Ca2, Cb  (Caton et al. 1982, mature H1 numbering)
# ANTIGENIC RIDGE:
#   H3: Koel et al. 2013 (Science 342:976-979)
#   H1: Xu et al. 2010 + Caton et al. 1982
#
# INPUT:   {aln_dir}/{gene}_aligned.fasta
# OUTPUT:  data/epitope/{gene}_{epitope,nonepitope,ridge,nonridge}_cols.txt
#          data/epitope/{gene}_site_map.csv
#
# FLAGS:
#   --gene X       run a single gene (default: H1 H3)
#   --aln-dir D    read alignments from D instead of PATHS$alignments
#                  (e.g. --aln-dir data/alignments/old while MAFFT is re-run)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
  library(Biostrings)
})

source(here("config.R"))

args    <- commandArgs(trailingOnly = TRUE)
genes   <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else c("H1", "H3")
aln_dir <- if ("--aln-dir" %in% args) args[which(args == "--aln-dir") + 1] else
           PATHS$alignments

out_dir <- PATHS$epitope
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

GENETIC <- as.list(Biostrings::GENETIC_CODE)

# ── Antigenic site definitions (mature protein numbering) ────────────────────
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
antigenic_sites_h1 <- list(
  Sa  = c(124:141),
  Sb  = c(156:168),
  Ca1 = c(166:177),
  Ca2 = c(203:211),
  Cb  = c(73:81)
)

# gene_meta: accession, signal peptide length, HA1/HA2 boundary, site sets, and
# a diagnostic motif from the mature N-terminus used to verify the fetched
# reference is actually the right subtype (see verify_reference_subtype).
gene_meta <- list(
  H3 = list(accession    = "EF614248.1",
            ref_name     = "ref_H3_Aichi68",
            signal_aa    = 16L,
            ha1_last_aa  = 329L,
            mature_motif = "PGNDNSTATLCLGHH",
            sites        = antigenic_sites_h3,
            ridge        = c(145, 155, 156, 158, 159, 189, 193)),
  H1 = list(accession    = "CY009444",
            ref_name     = "ref_H1_PR8",
            signal_aa    = 17L,
            ha1_last_aa  = 327L,
            mature_motif = "CIGYHANNST",
            sites        = antigenic_sites_h1,
            ridge        = c(133, 137, 138, 140, 141, 156, 157))
)

# ── Fetch reference (cache keyed by accession, so a changed accession refetches)
fetch_reference <- function(accession, ref_name, cache_dir) {
  tag        <- gsub("[^A-Za-z0-9]", "_", accession)
  cache_path <- file.path(cache_dir, paste0(ref_name, "_", tag, ".fasta"))
  if (file.exists(cache_path)) {
    message("  Using cached reference: ", basename(cache_path))
    return(cache_path)
  }
  message("  Fetching ", accession, " from NCBI GenBank...")
  seq <- tryCatch(read.GenBank(accession, as.character = TRUE),
                  error = function(e)
                    stop("Could not fetch ", accession, " from GenBank: ",
                         conditionMessage(e)))
  writeLines(c(paste0(">", ref_name), paste(seq[[1]], collapse = "")), cache_path)
  message("  Written: ", basename(cache_path))
  cache_path
}

# ── Translate a nucleotide vector in a given frame ───────────────────────────
translate_frame <- function(nt, offset) {
  starts <- seq(1 + offset, length(nt) - 2, by = 3)
  vapply(starts, function(s) {
    aa <- GENETIC[[paste0(nt[s], nt[s + 1], nt[s + 2])]]
    if (is.null(aa)) "?" else aa
  }, character(1))
}

# ── Locate the CDS: longest stop-free ORF from an ATG, across all frames ─────
# Scans all three frames rather than assuming one, and takes the longest ORF
# rather than the first ATG — both references carry 5' leaders and short
# spurious ORFs ahead of the real start codon.
locate_cds <- function(ref_nt) {
  frames  <- lapply(0:2, function(o) translate_frame(ref_nt, o))
  n_stops <- vapply(frames, function(a) sum(head(a, -1) == "*"), integer(1))
  best <- list(len = -1L)
  for (off in 0:2) {
    aa    <- frames[[off + 1L]]
    stops <- which(aa == "*")
    for (m in which(aa == "M")) {
      nxt <- stops[stops > m][1]
      end <- if (is.na(nxt)) length(aa) else nxt - 1L
      if (end - m + 1L > best$len)
        best <- list(len = end - m + 1L, offset = off, start = m, end = end,
                     aa = aa[m:end])
    }
  }
  if (best$len < 1L) stop("No ORF found in reference")
  message("    CDS: frame ", best$offset, ", ", best$len,
          " codons  (internal stops by frame: ", paste(n_stops, collapse = "/"), ")")
  if (best$len < 500L || best$len > 600L)
    warning("Reference ORF is ", best$len, " codons — expected ~560-570 for HA.")
  best
}

# ── Verify the fetched reference is the expected subtype ─────────────────────
# The original H3 accession returned an H1 HA, and nothing caught it because
# nothing ever translated the reference and looked at the protein. This checks
# the mature N-terminus against a subtype-diagnostic motif:
#   H1 mature: DT(I/L)CIGYHANNSTDTV...
#   H3 mature: Q(D/K)LPGNDNSTATLCLGHHAVPN...
verify_reference_subtype <- function(cds_aa, signal_aa, motif, gene) {
  mature <- paste(cds_aa[(signal_aa + 1):min(signal_aa + 40, length(cds_aa))],
                  collapse = "")
  if (!grepl(motif, mature, fixed = TRUE)) {
    stop("Reference for ", gene, " does not look like ", gene, ".\n",
         "  expected motif : ", motif, "\n",
         "  mature N-term  : ", mature, "\n",
         "  This is how the CY121680 mix-up (an H1 sequence used as the H3 ",
         "reference) went unnoticed. Check the accession.")
  }
  message("    subtype check OK (mature N-term: ", substr(mature, 1, 24), "...)")
}

# ── Add reference to the alignment and capture MAFFT's base→column map ───────
# --keeplength is deliberately NOT used (it deletes non-fitting reference
# bases). --mapout emits an explicit map to the ORIGINAL alignment's columns.
add_reference <- function(aln_path, ref_fasta, mafft_bin, threads = 1L) {
  tmp_out  <- tempfile(fileext = ".fasta")
  map_file <- paste0(ref_fasta, ".map")
  if (file.exists(map_file)) unlink(map_file)
  ret <- system2(mafft_bin,
                 args   = c("--quiet", "--thread", threads,
                            "--add", ref_fasta, "--mapout", aln_path),
                 stdout = tmp_out, stderr = FALSE)
  if (ret != 0 || !file.exists(tmp_out))
    stop("mafft --add failed (return code ", ret, ")")
  if (!file.exists(map_file))
    stop("mafft --mapout produced no map file: ", map_file)
  list(aln = tmp_out, map = map_file)
}

# ── Parse a mafft --mapout file ──────────────────────────────────────────────
# Format:  letter, position in the original sequence, position in the alignment
#          a, 1, 1
#          t, 7, -     <- insertion relative to the alignment; unmappable
# Returns an integer vector indexed by reference position -> alignment column
# (NA where unmappable).
read_mafft_map <- function(map_file) {
  lines <- readLines(map_file, warn = FALSE)
  lines <- lines[!grepl("^\\s*(>|#)", lines) & nzchar(trimws(lines))]
  parts <- strsplit(lines, ",")
  parts <- parts[lengths(parts) >= 3]
  if (!length(parts)) stop("No rows parsed from ", map_file)
  ref_pos <- as.integer(trimws(vapply(parts, `[`, character(1), 2)))
  aln_pos <- trimws(vapply(parts, `[`, character(1), 3))
  out <- rep(NA_integer_, max(ref_pos, na.rm = TRUE))
  ok  <- aln_pos != "-" & !is.na(ref_pos)
  out[ref_pos[ok]] <- as.integer(aln_pos[ok])
  out
}

# ── Build codon → alignment-column map ───────────────────────────────────────
# Rownames are mature aa positions (negative through the signal peptide).
# Numbering is assigned BEFORE unmappable codons are dropped, so residue
# positions stay anchored to the reference regardless of what is skipped.
build_codon_map <- function(ref_fasta, map_file, signal_aa, motif, gene) {
  ref_nt <- toupper(as.vector(read.dna(ref_fasta, format = "fasta",
                                       as.character = TRUE, as.matrix = TRUE)[1, ]))
  cds <- locate_cds(ref_nt)
  verify_reference_subtype(cds$aa, signal_aa, motif, gene)

  codon_pos <- t(vapply(cds$start:cds$end, function(k) {
    s <- 1 + cds$offset + (k - 1) * 3
    c(s, s + 1, s + 2)
  }, numeric(3)))

  orig2col <- read_mafft_map(map_file)
  n_unmapped <- sum(is.na(orig2col[seq_len(min(length(orig2col), length(ref_nt)))]))
  message("    mapout: ", n_unmapped, " of ", length(ref_nt),
          " reference bases have no alignment column")

  idx  <- ifelse(codon_pos <= length(orig2col), codon_pos, NA_integer_)
  cols <- matrix(orig2col[idx], ncol = 3)
  ok   <- stats::complete.cases(cols)
  if (!any(ok)) stop("No reference codons could be mapped for ", gene)
  if (sum(!ok) > 0)
    message("    ", sum(!ok), " of ", nrow(cols), " CDS codons unmappable, skipped")

  rownames(cols) <- as.character(seq_len(nrow(cols)) - signal_aa)
  cols <- cols[ok, , drop = FALSE]
  message("    mapped ", nrow(cols), " codons to alignment columns")
  cols
}

# ── Validate the codon map against real sequences ────────────────────────────
# Translates a sample of real sequences through the map and checks the result
# looks like a protein. Nothing in the original version ever did this, which is
# why both the wrong-subtype and frame bugs went unnoticed.
validate_codon_map <- function(aln, codon_map, signal_aa, gene,
                               n_sample = 50, max_stop_frac = 0.02) {
  set.seed(1)
  samp <- sample(rownames(aln), min(n_sample, nrow(aln)))
  start_key <- as.character(1L - signal_aa)

  stop_fracs <- numeric(0); met_start <- 0L
  for (s in samp) {
    aa <- apply(codon_map, 1, function(cc) {
      a <- GENETIC[[toupper(paste(aln[s, cc], collapse = ""))]]
      if (is.null(a)) NA_character_ else a
    })
    aa <- aa[!is.na(aa)]
    if (!length(aa)) next
    stop_fracs <- c(stop_fracs, mean(head(aa, -1) == "*"))
    # Check the codon at precursor residue 1, not merely the first retained
    # row: if codon 1 was skipped, row 1 is a later residue and a naive check
    # false-alarms.
    if (start_key %in% names(aa) && aa[[start_key]] == "M")
      met_start <- met_start + 1L
  }
  med <- stats::median(stop_fracs)
  message("    validation: median internal-stop fraction = ", signif(med, 3),
          "; ", met_start, "/", length(samp), " start with Met")
  if (med > max_stop_frac)
    stop("Codon map failed validation for ", gene, ": median internal-stop ",
         "fraction ", signif(med, 3), " exceeds ", max_stop_frac,
         " — the reading frame is probably mis-registered.")
  invisible(med)
}

# ── Resolve site positions to alignment columns ──────────────────────────────
sites_to_cols <- function(sites_list, codon_map) {
  bind_rows(lapply(names(sites_list), function(site_name) {
    positions <- sites_list[[site_name]]
    present   <- positions[as.character(positions) %in% rownames(codon_map)]
    if (!length(present)) return(NULL)
    data.frame(
      site          = site_name,
      residue_pos   = rep(present, each = 3),
      alignment_col = as.vector(t(codon_map[as.character(present), , drop = FALSE]))
    )
  }))
}

# ── Main loop ────────────────────────────────────────────────────────────────
for (gene in genes) {
  message("\n[", gene, "] Mapping epitope sites...")
  meta     <- gene_meta[[gene]]
  aln_path <- file.path(aln_dir, paste0(gene, "_aligned.fasta"))
  if (!file.exists(aln_path)) stop("Alignment not found: ", aln_path)

  ref_fasta <- fetch_reference(meta$accession, meta$ref_name, out_dir)

  message("  Running mafft --add --mapout (temp file)...")
  added <- add_reference(aln_path, ref_fasta, TOOLS$mafft,
                         threads = PARAMS$mafft_threads)

  # The map reports ORIGINAL-alignment columns, so only the original alignment
  # is needed from here — the extended one is not used.
  orig_aln    <- read.dna(aln_path, format = "fasta",
                          as.character = TRUE, as.matrix = TRUE)
  n_orig_cols <- ncol(orig_aln)
  message("  Alignment columns: ", n_orig_cols)

  codon_map <- build_codon_map(ref_fasta, added$map, meta$signal_aa,
                               meta$mature_motif, gene)
  validate_codon_map(orig_aln, codon_map, meta$signal_aa, gene)
  unlink(c(added$aln, added$map))

  # ── Epitope / ridge columns ────────────────────────────────────────────────
  epitope_df <- sites_to_cols(meta$sites, codon_map) |>
    filter(alignment_col <= n_orig_cols)
  epitope_cols    <- sort(unique(epitope_df$alignment_col))
  nonepitope_cols <- setdiff(seq_len(n_orig_cols), epitope_cols)

  ridge_present <- meta$ridge[as.character(meta$ridge) %in% rownames(codon_map)]
  ridge_cols    <- sort(unique(as.vector(
    codon_map[as.character(ridge_present), , drop = FALSE])))
  ridge_cols    <- ridge_cols[ridge_cols <= n_orig_cols]
  nonridge_cols <- setdiff(seq_len(n_orig_cols), ridge_cols)

  n_sites_def <- length(unique(unlist(meta$sites)))
  n_sites_map <- n_distinct(epitope_df$residue_pos)
  message("  Epitope: ", n_sites_map, "/", n_sites_def, " defined positions mapped",
          " (", length(epitope_cols), " columns)")
  message("  Ridge:   ", length(ridge_present), "/", length(meta$ridge),
          " positions mapped (", length(ridge_cols), " columns)")
  if (n_sites_map < n_sites_def)
    warning(gene, ": ", n_sites_def - n_sites_map,
            " defined antigenic positions could not be mapped")

  # ── Site map table ─────────────────────────────────────────────────────────
  site_map <- bind_rows(lapply(seq_len(nrow(codon_map)), function(i) {
    mature <- as.integer(rownames(codon_map)[i])
    cols   <- codon_map[i, ]
    cols   <- cols[cols <= n_orig_cols]
    if (!length(cols)) return(NULL)
    region <- if (mature <= 0) "signal" else
              if (mature <= meta$ha1_last_aa) "HA1" else "HA2"
    hit  <- names(Filter(function(s) mature %in% s, meta$sites))
    site <- if (length(hit)) hit[1] else NA_character_
    if (!is.na(site)) region <- paste0("epitope_", site)
    if (mature %in% meta$ridge) region <- paste0(region, "_ridge")
    data.frame(alignment_col = cols, nt_start = cols, nt_end = cols,
               residue_pos = mature, region = region, site = site,
               stringsAsFactors = FALSE)
  }))

  writeLines(as.character(epitope_cols),
             file.path(out_dir, paste0(gene, "_epitope_cols.txt")))
  writeLines(as.character(nonepitope_cols),
             file.path(out_dir, paste0(gene, "_nonepitope_cols.txt")))
  writeLines(as.character(ridge_cols),
             file.path(out_dir, paste0(gene, "_ridge_cols.txt")))
  writeLines(as.character(nonridge_cols),
             file.path(out_dir, paste0(gene, "_nonridge_cols.txt")))
  write_csv(site_map, file.path(out_dir, paste0(gene, "_site_map.csv")))
  message("  Written: ", gene, "_{epitope,nonepitope,ridge,nonridge}_cols.txt",
          " / _site_map.csv")
}

message("\nDone.")
