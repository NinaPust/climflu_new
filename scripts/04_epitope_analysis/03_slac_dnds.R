# scripts/04_epitope_analysis/03_slac_dnds.R
#
# PURPOSE: Per-terminal-branch dN and dS from HyPhy SLAC (Single-Likelihood
#          Ancestor Counting), for H1 and H3. Replaces the parsimony counting in
#          09_substitutions.R as the source for figure 1D.
#
# WHY SLAC:
#   - Maximum likelihood, not parsimony. SLAC fits a codon substitution model
#     (MG94xREV), reconstructs ancestral states under it, then counts. Parsimony
#     counts the minimum number of changes and cannot correct for multiple hits
#     or unequal mutation opportunity.
#   - Substitution counts arrive ALREADY NORMALISED by expected site counts
#     (ES/EN per codon site, from the fitted codon model). That is Reviewer 1's
#     normalisation comment answered directly, and more principled than the
#     Nei-Gojobori counts 09_substitutions.R uses, which ignore transition/
#     transversion bias and codon usage.
#   - No tree mapping. This is the step that killed the original yn00 route:
#     05_map_to_tree_nnls.R had to fit dN and dS branch lengths onto the
#     topology from a pairwise matrix, and failed. SLAC reconstructs ancestors
#     on the given tree directly, so nothing has to be mapped back.
#   - Reviewer 2 suggested it.
#
# ONE RUN PER GENE COVERS EVERY REGION:
#   SLAC's JSON stores substitution counts as per-branch x per-codon-site
#   vectors, and expected site counts (ES/EN) per codon site. So the whole CDS
#   is analysed once and regions are obtained afterwards by summing over the
#   relevant codons. Running SLAC separately per region would refit the codon
#   model on a subalignment each time and give inconsistent site counts.
#
# REGIONS (derived from the site map, so they inherit the reading-frame fix):
#   full     every CDS codon (signal + HA1 + HA2)
#   HA1      HA1 excluding antigenic sites (site map labels those epitope_*)
#   HA2      HA2 subunit — conserved, the internal negative control
#   epitope  all antigenic-site codons
#   ridge    the 7 antigenic-ridge codons (a subset of epitope)
#
# ZEROES: terminal branches with no synonymous change give dS = 0 and therefore
#   undefined omega; branches with no change at all give dN = dS = 0. SLAC counts
#   more accurately than parsimony but cannot create substitutions that are not
#   there, so the zero fraction stays high in the short regions. Counts, expected
#   sites, dN and dS are all written out unconditionally; how zeroes are handled
#   is a modelling decision left to the downstream script.
#
# INPUT:
#   data/epitope/{gene}_full_aligned.fasta   in-frame CDS from 02
#   data/epitope/{gene}_site_map.csv         codon -> region
#   data/trees/timed/{gene}_timetree.nwk               topology
#   data/rates/{gene}_rates_epitope.csv      to join onto (optional)
#
# OUTPUT:
#   data/slac/{gene}/aln.fasta, tree.nwk, slac.json, slac.log
#   data/slac/{gene}/stop_codons.csv     which sequences carried an
#                                                  internal stop, and whether it
#                                                  was masked or dropped
#   data/rates/{gene}_slac_branches.csv   all branches, all regions
#   data/rates/{gene}_rates_slac.csv      tips only, joined to rates
#
# FLAGS:
#   --gene X        run a single gene (default: H1 H3)
#   --threads N     HyPhy CPU cap (default 4; shared machine)
#   --samples N     ancestral-uncertainty resamples (default 10; HyPhy's own
#                   default is 100, which is slower and only widens the CIs,
#                   which nothing downstream reads)
#   --force         re-run SLAC even if slac.json already exists
#   --parse-only    skip HyPhy, parse existing slac.json
#   --estimator X   RESOLVED (default) or AVERAGED, for the whole-CDS by-branch
#                   cross-check table only
#   --max-stops N   sequences with more than N internal stop codons are dropped
#                   rather than masked (default 5). See mask_stop_codons().
#
# RUNTIME: measured linear in tip count at 4 cores — 150 tips 22s, 400 tips 59s,
#   800 tips 116s. Extrapolates to roughly 26 min for H1 (~10.9k tips) and
#   31 min for H3 (~13k tips). Run it in tmux.
#
# MEMORY: the JSON is large (~12 MB per 800 tips, so ~170 MB for H1 at full
#   size) and jsonlite expands it several-fold while parsing. Expect a peak of
#   1-2 GB per gene. The genes are processed one at a time for that reason; do
#   not run two copies of this script concurrently.

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(Biostrings)
  library(jsonlite)
  library(dplyr)
  library(readr)
})

source(here("config.R"))

args       <- commandArgs(trailingOnly = TRUE)
genes      <- if ("--gene"      %in% args) args[which(args == "--gene")      + 1] else c("H1", "H3")
n_thread   <- if ("--threads"   %in% args) as.integer(args[which(args == "--threads") + 1]) else 4L
n_samples  <- if ("--samples"   %in% args) as.integer(args[which(args == "--samples") + 1]) else 10L
force_run  <- "--force"      %in% args
parse_only <- "--parse-only" %in% args
estimator  <- if ("--estimator" %in% args) args[which(args == "--estimator") + 1] else "RESOLVED"
max_stops  <- if ("--max-stops" %in% args) as.integer(args[which(args == "--max-stops") + 1]) else 5L

stopifnot(estimator %in% c("RESOLVED", "AVERAGED"))

hyphy_bin <- if (!is.null(TOOLS$hyphy) && nzchar(TOOLS$hyphy)) TOOLS$hyphy else Sys.which("hyphy")
if (!nzchar(hyphy_bin) || !file.exists(hyphy_bin))
  stop("HyPhy not found. Set TOOLS$hyphy in config.R (expected e.g. ",
       "/home/nina/miniforge3/envs/hyphy_env/bin/hyphy)")

slac_root <- here("data", "slac")
dir.create(slac_root, recursive = TRUE, showWarnings = FALSE)


#' Map each CDS codon to a region label set.
#'
#' The site map is per nucleotide column. The in-frame CDS alignment written by
#' 02_split_alignments.R uses those columns in sorted order, so codon k is
#' columns (3k-2):(3k) of that sorted vector. A codon's region is taken from its
#' first column; all three columns of a codon share a region by construction.
#'
#' @param site_map_path Path to {gene}_site_map.csv.
#' @param n_codons      Number of codons in the CDS alignment (for validation).
#' @return list of integer vectors, one per region, holding codon indices.
codon_regions <- function(site_map_path, n_codons) {
  sm <- read_csv(site_map_path, show_col_types = FALSE)
  cols <- sort(unique(sm$alignment_col))
  if (length(cols) %% 3 != 0)
    stop("site map column count (", length(cols), ") is not a multiple of 3")
  if (length(cols) / 3 != n_codons)
    stop("site map implies ", length(cols) / 3, " codons but the alignment has ",
         n_codons, " — regenerate the alignment with 02_split_alignments.R")

  first_col <- cols[seq(1, length(cols), by = 3)]
  region    <- sm$region[match(first_col, sm$alignment_col)]

  out <- list(
    full    = seq_len(n_codons),
    HA1     = which(region == "HA1"),
    HA2     = which(grepl("^HA2", region)),
    epitope = which(grepl("^epitope", region)),
    ridge   = which(grepl("_ridge$", region))
  )
  out[lengths(out) > 0]
}


#' Mask internal stop codons, dropping sequences that carry too many.
#'
#' HyPhy refuses to start if the alignment contains any stop codon
#' ("ASSERTION FAILED: The input alignment ... must not contain stop codons").
#' A handful of GISAID sequences carry one, almost always a sequencing error
#' rather than a real read-through: in this dataset H1 has none and H3 has
#' exactly one stop codon in one sequence, out of 8.47 million codons.
#'
#' Masking the offending codon as NNN is preferred to dropping the sequence,
#' because every tip carries its own location, date and climate values and is a
#' data point in the downstream regression — losing a tip costs more than losing
#' one codon in one sequence. HyPhy reads NNN as missing and simply contributes
#' nothing at that position.
#'
#' A sequence with many stops is a different matter: that suggests a frameshift
#' or a pseudogene rather than a point error, and masking would quietly keep bad
#' data. Those are dropped and named in the log.
#'
#' @param aln       DNAStringSet, in-frame CDS.
#' @param max_stops Sequences with more than this many internal stops are
#'                  dropped rather than masked.
#' @param log_path  Where to write the per-sequence record (NULL to skip).
#' @return The cleaned DNAStringSet.
mask_stop_codons <- function(aln, max_stops = 5L, log_path = NULL) {
  gc_tab <- Biostrings::GENETIC_CODE
  m      <- as.matrix(aln)
  n_cod  <- ncol(m) / 3L

  codons_of <- function(i) apply(matrix(toupper(m[i, ]), nrow = 3L), 2L,
                                 paste0, collapse = "")

  stop_idx <- lapply(seq_len(nrow(m)), function(i) {
    aa <- gc_tab[codons_of(i)]
    which(!is.na(aa) & aa == "*")
  })
  n_stops <- lengths(stop_idx)

  if (all(n_stops == 0L)) {
    message("  stop codons: none")
    return(aln)
  }

  drop_i <- which(n_stops > max_stops)
  mask_i <- which(n_stops > 0L & n_stops <= max_stops)

  message("  stop codons: ", sum(n_stops), " in ", sum(n_stops > 0), " sequence(s)")
  if (length(mask_i))
    message("    masking ", sum(n_stops[mask_i]), " codon(s) as NNN in ",
            length(mask_i), " sequence(s)")
  if (length(drop_i))
    message("    DROPPING ", length(drop_i), " sequence(s) with > ", max_stops,
            " stops: ", paste(head(rownames(m)[drop_i], 5), collapse = ", "),
            if (length(drop_i) > 5) " ..." else "")

  if (!is.null(log_path)) {
    rec <- do.call(rbind, lapply(which(n_stops > 0L), function(i) {
      data.frame(sequence = rownames(m)[i], n_stops = n_stops[i],
                 codon_index = paste(stop_idx[[i]], collapse = ";"),
                 action = if (i %in% drop_i) "dropped" else "masked",
                 stringsAsFactors = FALSE)
    }))
    utils::write.csv(rec, log_path, row.names = FALSE)
    message("    record written: ", basename(log_path))
  }

  for (i in mask_i) {
    for (k in stop_idx[[i]]) m[i, ((k - 1L) * 3L + 1L):(k * 3L)] <- "N"
  }

  keep <- setdiff(seq_len(nrow(m)), drop_i)
  out  <- DNAStringSet(apply(m[keep, , drop = FALSE], 1L, paste0, collapse = ""))
  names(out) <- rownames(m)[keep]
  out
}


#' Prepare the SLAC input files for one gene and return their paths.
prepare_inputs <- function(gene, gdir) {
  aln_path <- file.path(PATHS$epitope, paste0(gene, "_full_aligned.fasta"))
  tree_in  <- file.path(here(), "data", "trees", "timed",
                        paste0(gene, "_timetree.nwk"))
  for (p in c(aln_path, tree_in))
    if (!file.exists(p)) stop("Required file not found: ", p)

  aln <- readDNAStringSet(aln_path)
  tr  <- read.tree(tree_in)

  common <- intersect(names(aln), tr$tip.label)
  if (length(common) < 50)
    stop(gene, ": only ", length(common), " tips shared between alignment and tree")
  message("  alignment: ", length(aln), " seqs  tree: ", Ntip(tr),
          " tips  shared: ", length(common))

  aln <- aln[common]

  # HyPhy aborts on any stop codon, so clean them before writing the input.
  aln <- mask_stop_codons(aln, max_stops = max_stops,
                          log_path = file.path(gdir, "stop_codons.csv"))
  if (length(aln) < length(common)) {
    common <- names(aln)
    message("  tips after stop-codon filtering: ", length(common))
  }

  tr  <- keep.tip(tr, common)
  # HyPhy treats {...} in node labels as branch annotation; strip them so the
  # tree is read as a plain topology.
  tr$node.label <- NULL

  f_aln  <- file.path(gdir, "aln.fasta")
  f_tree <- file.path(gdir, "tree.nwk")
  writeXStringSet(aln, f_aln)
  write.tree(tr, f_tree)
  list(aln = f_aln, tree = f_tree, n_tips = length(common),
       n_codons = unique(width(aln))[1] / 3)
}


#' Parse a SLAC JSON into per-branch, per-region dN/dS.
#'
#' Counts come from the per-branch x per-codon substitution vectors in
#' `branch attributes`.
#'
#' EXPECTED SITE COUNTS ARE PER BRANCH, NOT CONSTANT. They depend on how much of
#' the sequence is scoreable on that branch: a branch whose sequences carry gaps
#' or Ns contributes no expected sites at those codons. In H1 this is a large
#' effect — ES ranges 364.7 to 527.6 and EN 782.7 to 1123.8 across branches,
#' driven almost entirely by missing data (e.g. EPI_ISL_66521 has 407 of 565
#' codons scoreable and correspondingly ~30% fewer expected sites).
#'
#' Using a single constant per region, as an earlier version of this function
#' did, divides gappy branches by a denominator that is too large and so
#' understates their dN and dS by up to ~30% — a bias that tracks data quality,
#' which in turn tracks collection era and region.
#'
#' SLAC publishes per-branch ES/EN only for the WHOLE alignment
#' (`MLE$content$<part>$by-branch`). For sub-regions this function reconstructs
#' them as the sum of per-codon ES/EN over the region's codons, restricted to
#' codons that are scoreable on that branch. Validated against SLAC's own
#' whole-gene values over 2000 H1 branches: median difference +0.29 on ES
#' (base ~524) and -2.05 on EN (base ~1119), i.e. under 0.2%. The residual is
#' because SLAC averages the site counts of the two codons at each end of a
#' branch, whereas this uses the branch's own reconstructed codon.
#'
#' Both are reported for `full`: the reconstruction (`EN`/`ES`, consistent with
#' the sub-regions) and HyPhy's authoritative value (`EN_hyphy`/`ES_hyphy`), so
#' the agreement can be checked directly in the output.
#'
#' @param json_path Path to slac.json.
#' @param regions   Output of codon_regions().
#' @param estimator "RESOLVED" or "AVERAGED" (which ES/EN tables to use).
#' @return data frame, one row per branch per region.
parse_slac <- function(json_path, regions, estimator = "RESOLVED") {
  message("  parsing JSON (", round(file.size(json_path) / 1048576, 1), " MB)...")
  j <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)

  part <- names(j$`branch attributes`)
  part <- part[part != "attributes"][1]
  ba   <- j$`branch attributes`[[part]]

  hdr  <- vapply(j$MLE$headers, `[`, character(1), 1)
  ES_i <- match("ES", hdr); EN_i <- match("EN", hdr)
  if (is.na(ES_i) || is.na(EN_i)) stop("ES/EN columns not found in MLE headers")

  bysite <- j$MLE$content[[part]]$`by-site`[[estimator]]
  if (is.null(bysite)) stop("by-site table '", estimator, "' not found in JSON")
  ES_site <- as.numeric(bysite[, ES_i])
  EN_site <- as.numeric(bysite[, EN_i])

  # HyPhy's authoritative per-branch whole-gene ES/EN, keyed by branch name.
  bybr <- j$MLE$content[[part]]$`by-branch`
  hy_ES <- hy_EN <- NULL
  if (!is.null(bybr) && !is.null(bybr[[estimator]])) {
    bnames <- as.character(unlist(bybr$NAMES))
    btab   <- bybr[[estimator]]
    hy_ES  <- setNames(as.numeric(btab[, ES_i]), bnames)
    hy_EN  <- setNames(as.numeric(btab[, EN_i]), bnames)
    message("  by-branch ES/EN available for ", length(bnames), " branches")
  } else {
    warning("by-branch table absent; whole-gene cross-check unavailable")
  }

  node_names <- names(ba)
  message("  branches with attributes: ", length(node_names))

  rows <- lapply(node_names, function(nm) {
    v  <- ba[[nm]]
    ns <- v[["nonsynonymous substitution count"]]
    sy <- v[["synonymous substitution count"]]
    if (is.null(ns) || is.null(sy)) return(NULL)
    # Stored as a length-1 list holding the per-codon vector.
    ns <- as.numeric(unlist(ns)); sy <- as.numeric(unlist(sy))
    if (length(ns) != length(ES_site)) return(NULL)

    # Scoreable codons on this branch: no gap, no ambiguity.
    cod <- unlist(v[["codon"]])
    scoreable <- if (!is.null(cod) && length(cod) == length(ES_site)) {
      !grepl("[-?]", cod) & !grepl("N", cod, ignore.case = TRUE)
    } else {
      rep(TRUE, length(ES_site))   # no codon labels: fall back to all sites
    }

    orig  <- v[["original name"]]
    label <- if (!is.null(orig) && nzchar(orig)) orig else nm

    do.call(rbind, lapply(names(regions), function(rg) {
      idx <- regions[[rg]]
      idx <- idx[scoreable[idx]]
      data.frame(
        node     = label,
        region   = rg,
        N        = sum(ns[regions[[rg]]]),
        S        = sum(sy[regions[[rg]]]),
        EN       = sum(EN_site[idx]),
        ES       = sum(ES_site[idx]),
        n_codons_scored = length(idx),
        EN_hyphy = if (rg == "full" && !is.null(hy_EN)) unname(hy_EN[nm]) else NA_real_,
        ES_hyphy = if (rg == "full" && !is.null(hy_ES)) unname(hy_ES[nm]) else NA_real_,
        stringsAsFactors = FALSE
      )
    }))
  })

  out <- bind_rows(rows[!vapply(rows, is.null, logical(1))])
  if (!nrow(out)) stop("No branch rows parsed from ", json_path)

  # Validation: how closely does the reconstruction match HyPhy's own values?
  chk <- out[out$region == "full" & !is.na(out$EN_hyphy), , drop = FALSE]
  if (nrow(chk) > 0) {
    message("  whole-gene ES/EN reconstruction vs HyPhy (median diff): ",
            "ES ", sprintf("%+.2f", median(chk$ES - chk$ES_hyphy)),
            "  EN ", sprintf("%+.2f", median(chk$EN - chk$EN_hyphy)),
            "   (on bases of ~", round(median(chk$ES_hyphy)), " / ~",
            round(median(chk$EN_hyphy)), ")")
  }

  out |>
    mutate(
      dN    = ifelse(EN > 0, N / EN, NA_real_),
      dS    = ifelse(ES > 0, S / ES, NA_real_),
      omega = ifelse(!is.na(dS) & dS > 0, dN / dS, NA_real_)
    )
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for (gene in genes) {
  message("\n========== ", gene, " ==========")
  gdir <- file.path(slac_root, gene)
  dir.create(gdir, recursive = TRUE, showWarnings = FALSE)
  json_path <- file.path(gdir, "slac.json")

  prep <- prepare_inputs(gene, gdir)
  regions <- codon_regions(
    file.path(PATHS$epitope, paste0(gene, "_site_map.csv")), prep$n_codons)
  message("  regions (codons): ",
          paste(sprintf("%s=%d", names(regions), lengths(regions)), collapse = "  "))

  # ── Run SLAC ────────────────────────────────────────────────────────────────
  if (parse_only) {
    if (!file.exists(json_path)) { warning(gene, ": no slac.json to parse"); next }
    message("  --parse-only: using existing slac.json")
  } else if (file.exists(json_path) && !force_run) {
    message("  slac.json exists, skipping HyPhy (use --force to re-run)")
  } else {
    if (file.exists(json_path)) unlink(json_path)
    message("  running SLAC on ", prep$n_tips, " tips, ", prep$n_codons,
            " codons, ", n_thread, " threads...")
    message("  (expect roughly ", round(prep$n_tips * 0.145 / 60), " min; ",
            "progress goes to ", file.path(gdir, "slac.log"), ")")
    t0  <- proc.time()
    ret <- system2(
      hyphy_bin,
      args = c(paste0("CPU=", n_thread), "slac",
               "--alignment", shQuote(prep$aln),
               "--tree",      shQuote(prep$tree),
               "--samples",   n_samples,
               "--output",    shQuote(json_path)),
      stdout = file.path(gdir, "slac.log"), stderr = file.path(gdir, "slac.log")
    )
    el <- round((proc.time() - t0)[["elapsed"]] / 60, 1)
    if (ret != 0 || !file.exists(json_path)) {
      warning(gene, ": SLAC failed (exit ", ret, ") — see ",
              file.path(gdir, "slac.log"))
      next
    }
    message("  SLAC finished in ", el, " min")
  }

  # ── Parse ───────────────────────────────────────────────────────────────────
  br <- parse_slac(json_path, regions, estimator)
  br$gene <- gene
  write_csv(br, file.path(PATHS$rates, paste0(gene, "_slac_branches.csv")))
  message("  Written: ", gene, "_slac_branches.csv  (", nrow(br), " rows)")

  # ── Tip-level wide table, joined onto the rates CSV ─────────────────────────
  tip_labels <- read.tree(prep$tree)$tip.label
  brt <- br[br$node %in% tip_labels, , drop = FALSE]
  message("  terminal branches matched: ", n_distinct(brt$node), " of ",
          length(tip_labels))

  wide <- brt |>
    select(node, region, N, S, EN, ES, n_codons_scored, dN, dS, omega) |>
    tidyr::pivot_wider(
      names_from  = region,
      values_from = c(N, S, EN, ES, n_codons_scored, dN, dS, omega),
      names_glue  = "slac_{region}_{.value}"
    ) |>
    rename(accession = node)

  for (rg in names(regions)) {
    for (q in c("dN", "dS", "omega")) {
      src <- paste0("slac_", rg, "_", q)
      if (src %in% names(wide))
        wide[[paste0("log_", src)]] <-
          ifelse(!is.na(wide[[src]]) & wide[[src]] > 0, log(wide[[src]]), NA_real_)
    }
  }

  base_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))
  if (!file.exists(base_path))
    base_path <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))
  if (file.exists(base_path)) {
    base <- read_csv(base_path, show_col_types = FALSE)
    base <- base |> select(-any_of(setdiff(names(wide), "accession")))
    out  <- left_join(base, wide, by = "accession")
  } else {
    warning(gene, ": no rates CSV to join onto; writing SLAC columns alone")
    out <- wide
  }
  out_path <- file.path(PATHS$rates, paste0(gene, "_rates_slac.csv"))
  write_csv(out, out_path)
  message("  Written: ", basename(out_path), "  (", nrow(out), " rows)")

  # ── Summary ─────────────────────────────────────────────────────────────────
  message("\n  Terminal-branch summary (", gene, "):")
  brt |>
    group_by(region) |>
    summarise(
      tips      = n(),
      `N>0`     = sprintf("%.1f%%", 100 * mean(N > 0)),
      `S>0`     = sprintf("%.1f%%", 100 * mean(S > 0)),
      `both>0`  = sprintf("%.1f%%", 100 * mean(N > 0 & S > 0)),
      `no change` = sprintf("%.1f%%", 100 * mean(N == 0 & S == 0)),
      mean_dN   = signif(mean(dN, na.rm = TRUE), 3),
      mean_dS   = signif(mean(dS, na.rm = TRUE), 3),
      `EN/ES`   = signif(mean(EN) / mean(ES), 3),
      `median EN` = signif(median(EN), 4),
      `min EN`    = signif(min(EN), 4),
      .groups = "drop"
    ) |>
    as.data.frame() |>
    print(row.names = FALSE)
}

message("\nDone.")
