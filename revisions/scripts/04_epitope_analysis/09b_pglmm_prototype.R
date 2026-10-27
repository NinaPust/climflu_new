# revisions/scripts/04_epitope_analysis/09b_pglmm_prototype.R
#
# PURPOSE: Prototype of the "Path 3" count-model approach on a subsample, to
#          check feasibility/runtime before committing to it for all genes and
#          regions.
#
#          Path 3 keeps the existing parsimony ancestral reconstruction and its
#          NS counts untouched (from 09_substitutions.R). It only changes what
#          happens to those counts afterwards:
#
#            current : log(ns_count / branch_time)        -> phylolm
#                      (zeros become NA and are DROPPED)
#            path 3  : ns_count as a Poisson response,
#                      log(N_sites x branch_time) as exposure -> pglmm
#                      (zeros are RETAINED as real observations)
#
#          Both target the same quantity: substitutions per site per year.
#          The difference is that path 3 keeps the ~88% of epitope tips that
#          currently get discarded.
#
# NOTE ON OFFSETS: phyr::pglmm does not support offset() terms (confirmed: no
#          offset handling anywhere in its source). Rather than silently losing
#          the exposure term, log_exposure is entered as an ordinary covariate.
#          If its coefficient is ~1 the offset assumption is supported; this is
#          a standard, slightly more flexible alternative to a fixed offset.
#
# CAVEAT:  Uses the existing epitope NS counts, which still contain the
#          gap-codon-scored-as-no-change issue. The epitope region is NOT
#          affected by the full/HA2 reading-frame bug (it is reference-anchored),
#          so these counts are usable for a feasibility test. Numbers here are
#          for checking the machinery, not for reporting.
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv       ({region}_ns_count column)
#   data/epitope/{gene}_{region}_aligned.fasta
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/phylolm/{gene}_{region}_pglmm_prototype.csv
#
# FLAGS:
#   --gene X     default H1
#   --region X   default epitope
#   --n N        number of tips to subsample (default 200; 0 = use all)
#   --full       use the full covariate set (default: scaled_temp only, since
#                200 tips with ~88% zeros cannot support 17 parameters)
#   --seed N     default 42

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
  library(phyr)
  library(phylolm)
  library(Biostrings)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "codon_sites.R"))

args   <- commandArgs(trailingOnly = TRUE)
gene   <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else "H1"
region <- if ("--region" %in% args) args[which(args == "--region") + 1] else "epitope"
n_sub  <- if ("--n"      %in% args) as.integer(args[which(args == "--n") + 1]) else 200L
seed   <- if ("--seed"   %in% args) as.integer(args[which(args == "--seed") + 1]) else 42L
use_full_model <- "--full" %in% args

set.seed(seed)
message("=== pglmm prototype: ", gene, " / ", region,
        "  (n = ", if (n_sub > 0) n_sub else "all", ") ===")

# ── Load ──────────────────────────────────────────────────────────────────────
rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))
aln_path   <- file.path(PATHS$epitope, paste0(gene, "_", region, "_aligned.fasta"))
tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
for (p in c(rates_path, aln_path, tree_path))
  if (!file.exists(p)) stop("Missing input: ", p)

df   <- read_csv(rates_path, show_col_types = FALSE)
tree <- read.tree(tree_path)
aln  <- read.dna(aln_path, format = "fasta", as.character = TRUE, as.matrix = TRUE)

count_col <- paste0(region, "_ns_count")
if (!count_col %in% names(df)) stop("Count column not found: ", count_col)

# ── Terminal branch times (from the FULL tree, before any subsetting) ─────────
# Pruning a tree changes terminal branch lengths, so exposure must be taken
# from the complete tree to stay comparable with the real analysis.
tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
tip_bl_full  <- setNames(tree$edge.length[tip_edge_idx], tree$tip.label)

# ── Per-tip expected non-synonymous site counts ───────────────────────────────
message("Computing Nei-Gojobori site counts...")
site_tab <- codon_site_table()
common   <- Reduce(intersect, list(tree$tip.label, rownames(aln), df$accession))
common   <- common[tip_bl_full[common] > 0]
message("  tips available (tree n aln n rates, positive branch length): ", length(common))

if (n_sub > 0 && n_sub < length(common)) common <- sample(common, n_sub)

sites <- t(vapply(common, function(a) seq_site_counts(aln[a, ], site_tab),
                  numeric(3)))
colnames(sites) <- c("S_sites", "N_sites", "n_codons_scored")

dat <- tibble(accession = common) |>
  left_join(df, by = "accession") |>
  mutate(
    ns_count     = .data[[count_col]],
    N_sites      = sites[accession, "N_sites"],
    branch_time  = tip_bl_full[accession],
    log_exposure = log(N_sites * branch_time)
  ) |>
  filter(!is.na(ns_count), is.finite(log_exposure), N_sites > 0)

message("  modelling n = ", nrow(dat),
        " | zeros = ", sum(dat$ns_count == 0),
        " (", round(100 * mean(dat$ns_count == 0), 1), "%)",
        " | nonzero = ", sum(dat$ns_count > 0))
message("  mean N_sites = ", round(mean(dat$N_sites), 1),
        " | mean branch_time = ", signif(mean(dat$branch_time), 3))

# ── Model specification ───────────────────────────────────────────────────────
covars <- if (use_full_model) {
  intersect(c("scaled_temp", "scaled_humidity", "scaled_delta_temp",
              "scaled_delta_humidity", "scaled_sampling_density",
              "scaled_pop_density"), names(dat))
} else {
  "scaled_temp"
}
message("  covariates: ", paste(covars, collapse = " + "))

dat <- dat |> filter(if_all(all_of(covars), is.finite))
sub_tree <- keep.tip(tree, dat$accession)
dat <- dat |> arrange(match(accession, sub_tree$tip.label))
dat$sp <- factor(dat$accession, levels = sub_tree$tip.label)

# ── Fit 1: phylogenetic Poisson count model (Path 3) ─────────────────────────
# (1|sp__) gives both an iid tip effect (overdispersion) and a phylogenetically
# structured tip effect. add.obs.re = FALSE avoids confounding the iid term with
# an additional observation-level effect, since there is one row per tip.
f_pglmm <- as.formula(paste("ns_count ~", paste(c(covars, "log_exposure"), collapse = " + "),
                            "+ (1|sp__)"))
message("\n[pglmm] ", deparse(f_pglmm))
t0 <- proc.time()
fit_pglmm <- tryCatch(
  pglmm(f_pglmm, data = dat, family = "poisson",
        cov_ranef = list(sp = sub_tree), add.obs.re = FALSE),
  error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL }
)
t_pglmm <- round((proc.time() - t0)[["elapsed"]], 1)
message("  elapsed: ", t_pglmm, "s")

# ── Fit 2: current approach (phylolm on log rate, nonzero subset only) ───────
dat_nz <- dat |> filter(ns_count > 0) |>
  mutate(log_ns_rate = log(ns_count / branch_time))
message("\n[phylolm, current approach] n = ", nrow(dat_nz),
        " (", nrow(dat) - nrow(dat_nz), " zero-count tips dropped)")
fit_phylolm <- NULL
if (nrow(dat_nz) >= 10) {
  tree_nz <- keep.tip(tree, dat_nz$accession)
  dat_nz  <- dat_nz |> arrange(match(accession, tree_nz$tip.label)) |> as.data.frame()
  rownames(dat_nz) <- dat_nz$accession
  f_pl <- as.formula(paste("log_ns_rate ~", paste(covars, collapse = " + ")))
  t1 <- proc.time()
  fit_phylolm <- tryCatch(
    phylolm(f_pl, data = dat_nz, phy = ape::unroot(tree_nz),
            model = PARAMS$phylolm_method),
    error = function(e) { message("  FAILED: ", conditionMessage(e)); NULL }
  )
  message("  elapsed: ", round((proc.time() - t1)[["elapsed"]], 1), "s")
} else {
  message("  too few nonzero tips to fit")
}

# ── Report ────────────────────────────────────────────────────────────────────
tidy_pglmm <- function(fit) {
  if (is.null(fit)) return(NULL)
  b <- as.data.frame(fit$B)
  data.frame(model = "pglmm_poisson", term = rownames(fit$B),
             estimate = b[[1]], std_error = fit$B.se[, 1],
             p.value = fit$B.pvalue[, 1], row.names = NULL)
}
tidy_phylolm <- function(fit) {
  if (is.null(fit)) return(NULL)
  cf <- as.data.frame(summary(fit)$coefficients)
  data.frame(model = "phylolm_lograte_nonzero", term = rownames(cf),
             estimate = cf[[1]], std_error = cf[[2]],
             p.value = cf[[4]], row.names = NULL)
}

out <- bind_rows(tidy_pglmm(fit_pglmm), tidy_phylolm(fit_phylolm))
if (!is.null(out) && nrow(out) > 0) {
  out$gene <- gene; out$region <- region
  out$n_total <- nrow(dat); out$n_zero <- sum(dat$ns_count == 0)
  cat("\n================ RESULTS ================\n")
  print(out |> mutate(across(c(estimate, std_error), \(x) round(x, 4)),
                      p.value = signif(p.value, 3)) |> as.data.frame(),
        row.names = FALSE)

  dir.create(PATHS$rates_phylolm, recursive = TRUE, showWarnings = FALSE)
  out_path <- file.path(PATHS$rates_phylolm,
                        paste0(gene, "_", region, "_pglmm_prototype.csv"))
  write_csv(out, out_path)
  message("\nWritten: ", out_path)
}

if (!is.null(fit_pglmm)) {
  cat("\n--- pglmm variance components ---\n")
  print(fit_pglmm$s2r)
  le <- out[out$model == "pglmm_poisson" & out$term == "log_exposure", "estimate"]
  if (length(le) == 1)
    cat("\nlog_exposure coefficient =", round(le, 3),
        "  (~1 supports treating it as an offset)\n")
}
message("\nDone. pglmm fit time: ", t_pglmm, "s for n = ", nrow(dat))
