# revisions/scripts/04_epitope_analysis/12a_phyloglm_any_change.R
#
# HURDLE PART 1 — "did any substitution occur on this terminal branch?"
#
# PURPOSE: Phylogenetic logistic regression of a binary change indicator against
#          climate predictors, for all four HA regions. This is the RATE half of
#          the hurdle: it asks whether change happened at all, not what kind.
#
#          Response: 1 if (ns_count + s_count) > 0, else 0.
#
#          Every tip gets an answer, so NOTHING IS DROPPED for being zero. This
#          is the direct answer to Reviewer 1's truncation-bias comment: the
#          submitted analysis modelled log(count / branch_time) and discarded
#          every tip with a zero numerator, which for the ridge was ~98% of them.
#
# EXPOSURE: a longer terminal branch, or one with more scoreable (non-gap)
#          codons, is mechanically more likely to carry a change whatever the
#          climate. So log((N_sites + S_sites) * branch_time) enters as a
#          covariate. Note it uses TOTAL sites, not just non-synonymous ones,
#          because the response counts synonymous changes too.
#
#          phyloglm has no offset() support, so exposure is an estimated term
#          rather than a fixed offset. Under a clean Poisson-like process its
#          coefficient should sit near 1; a value far from 1 is a finding about
#          the process, and 12c refits the same data as a count model to see
#          whether the deviation is an artefact of collapsing counts to binary.
#
#          Exposure is deliberately kept out of the FDR correction and the
#          forest plots: it is a measurement construct, not a hypothesis.
#
# INTERPRETATION: a positive temperature coefficient means change was MORE
#          LIKELY at higher temperatures, holding exposure and the other
#          predictors fixed. It says nothing about how much change or what kind
#          — a branch with one substitution and a branch with five are both
#          just "yes". That resolution is what Part 2 (12b) recovers.
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv   ({region}_ns_count, _s_count,
#                                          _N_sites, _S_sites from 09)
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/phylolm/{gene}_phyloglm_anychange_results.csv
#   data/rates/phylolm/all_phyloglm_anychange_results.csv
#   figures/epitope/phyloglm_anychange_forest.png
#   figures/epitope/phyloglm_anychange_forest_fdr.png
#
# FLAGS:
#   --gene X      run a single gene (default: H1 H3)
#   --region X    run a single region (default: full HA2 epitope ridge)
#   --no-year     drop the as.factor(year) term. Worth having for the ridge,
#                 where ~200-270 events against ~18 parameters is thin; the
#                 year factor alone costs 10 of them.
#   --fit-null    also fit a year-only model per region for delta-AIC. This
#                 DOUBLES runtime and nothing downstream currently reads it.
#   --threads N   BLAS/OMP cap (default 4; shared machine)
#   --max-n N     CALIBRATION ONLY — randomly subsample to N tips per fit.
#                 For checking the script runs and estimating runtime before
#                 committing to the full fit. Results from a subsample are NOT
#                 reportable: the filename is suffixed _subN so they cannot be
#                 confused with the real output.
#   --seed N      subsample seed (default 42)
#
# RUNTIME: phyloglm cost rises steeply with tip count (roughly n^2.5-3). All
#   four regions here use every tip (~9.3k H1, ~12.8k H3), so this is the
#   expensive script of the three. Run it in tmux.
#
#   MEMORY IS NOT A CONSTRAINT. phylolm/phyloglm use Ho & Ane's linear-time
#   algorithm and never materialise the dense n x n phylogenetic covariance, so
#   footprint stays flat: measured at ~290 MB for H1 with all tips. Running H1
#   and H3 concurrently is fine on memory; the only reason to serialise is the
#   4-thread CPU budget (two concurrent runs = ~8 cores).

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "phylolm_utils.R"))
source(here("revisions", "R", "hurdle_utils.R"))
source(here("revisions", "R", "forest_utils.R"))

args     <- commandArgs(trailingOnly = TRUE)
genes    <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
regions  <- if ("--region" %in% args) args[which(args == "--region") + 1] else
            c("full", "HA2", "epitope", "ridge")
use_year <- !("--no-year"  %in% args)
fit_null <- "--fit-null"   %in% args
n_thread <- if ("--threads" %in% args) as.integer(args[which(args == "--threads") + 1]) else 4L
max_n    <- if ("--max-n"   %in% args) as.integer(args[which(args == "--max-n")  + 1]) else NA_integer_
seed     <- if ("--seed"    %in% args) as.integer(args[which(args == "--seed")   + 1]) else 42L

cap_blas_threads(n_thread)

if (!is.na(max_n)) {
  set.seed(seed)
  message("\n*** CALIBRATION MODE: subsampling to ", max_n, " tips per fit. ***")
  message("*** Output is suffixed _sub", max_n, " and is NOT reportable. ***")
}
sfx <- if (is.na(max_n)) "" else paste0("_sub", max_n)

predictor_terms <- c(
  "scaled_temp", "scaled_humidity",
  "scaled_delta_temp", "scaled_delta_humidity",
  "scaled_sampling_density", "scaled_pop_density"
)

out_dir <- file.path(PATHS$rates_phylolm)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Fit per gene x region ─────────────────────────────────────────────────────

for (gene in genes) {
  message("\n[", gene, "] phyloglm — any-change (hurdle part 1)")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
  if (!file.exists(rates_path)) { warning("Rates not found: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Tree not found: ",  tree_path);  next }

  df        <- read_csv(rates_path, show_col_types = FALSE)
  full_tree <- read.tree(tree_path)

  # Branch durations from the UNPRUNED tree — pruning reassigns them.
  branch_time <- terminal_branch_times(full_tree)

  gene_results <- list()

  for (r in regions) {
    ns_col <- paste0(r, "_ns_count"); s_col  <- paste0(r, "_s_count")
    N_col  <- paste0(r, "_N_sites");  S_col  <- paste0(r, "_S_sites")
    needed <- c(ns_col, s_col, N_col, S_col)
    if (!all(needed %in% names(df))) {
      message("  [", r, "] missing columns (",
              paste(setdiff(needed, names(df)), collapse = ", "),
              ") — re-run 09_substitutions.R; skipping")
      next
    }

    sub <- df
    sub$any_change   <- as.integer(sub[[ns_col]] + sub[[s_col]] > 0)
    sub$log_exposure <- log((sub[[N_col]] + sub[[S_col]]) *
                            branch_time[sub$accession])

    keep <- c("accession", "any_change", "log_exposure", predictor_terms,
              if (use_year) "year")
    keep <- intersect(keep, names(sub))
    sub  <- sub[, keep, drop = FALSE]
    sub  <- sub[stats::complete.cases(sub) &
                is.finite(sub$log_exposure), , drop = FALSE]

    if (nrow(sub) < 50) { message("  [", r, "] <50 usable tips, skipping"); next }

    if (!is.na(max_n) && nrow(sub) > max_n)
      sub <- sub[sample(nrow(sub), max_n), , drop = FALSE]

    prepped  <- prep_tree_data(sub, full_tree)
    sub_df   <- prepped$data
    sub_tree <- prepped$tree

    n_events <- sum(sub_df$any_change == 1)
    message("  [", r, "] n = ", nrow(sub_df), "  events = ", n_events,
            " (", round(100 * n_events / nrow(sub_df), 1), "%)")
    if (n_events < 20 || n_events > nrow(sub_df) - 20) {
      message("    too few events (or non-events) for a stable fit, skipping")
      next
    }

    preds     <- intersect(predictor_terms, names(sub_df))
    year_term <- if (use_year && "year" %in% names(sub_df)) "as.factor(year)" else NULL
    n_par     <- length(preds) + 1L + (if (!is.null(year_term))
                                         length(unique(sub_df$year)) - 1L else 0L) + 1L
    if (n_events / n_par < 10)
      message("    NOTE: ", round(n_events / n_par, 1), " events per parameter ",
              "(<10) — estimates here will be unstable; consider --no-year")

    fml <- as.formula(paste("any_change ~",
                            paste(c(preds, "log_exposure", year_term),
                                  collapse = " + ")))

    t0  <- proc.time()
    res <- tryCatch(
      run_phyloglm(formula = fml, data = sub_df, tree = sub_tree,
                   subset = r, response = "any_change", gene = gene,
                   fit_null = fit_null),
      error = function(e) { message("    FAILED: ", conditionMessage(e)); NULL }
    )
    if (is.null(res)) next
    message("    fitted in ", round((proc.time() - t0)[["elapsed"]]), "s")

    exp_row <- res[res$term == "log_exposure", ]
    if (nrow(exp_row) == 1)
      message("    log_exposure coefficient = ", signif(exp_row$estimate, 3),
              " (Poisson expectation ~1; see 12c)")

    gene_results[[r]] <- res
  }

  if (length(gene_results) == 0) { message("  no successful fits"); next }
  gene_df <- bind_rows(gene_results)
  gene_out <- file.path(out_dir, paste0(gene, "_phyloglm_anychange_results", sfx, ".csv"))
  write_csv(gene_df, gene_out)
  message("  Written: ", basename(gene_out))
}

# ── Combine, FDR-correct, report ──────────────────────────────────────────────

all_res <- lapply(genes, function(g) {
  p <- file.path(out_dir, paste0(g, "_phyloglm_anychange_results", sfx, ".csv"))
  if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL
}) |> bind_rows()

if (nrow(all_res) == 0) stop("No results produced.")

all_res <- add_fdr(all_res, predictor_terms)
combined <- file.path(out_dir, paste0("all_phyloglm_anychange_results", sfx, ".csv"))
write_csv(all_res, combined)
message("\nWritten: ", basename(combined))

message("\nSignificant after BH-FDR (p.fdr < 0.05):")
sig <- all_res |>
  filter(term %in% predictor_terms, !is.na(p.fdr), p.fdr < 0.05) |>
  select(gene, subset, term, estimate, p.value, p.fdr, n, n_events) |>
  mutate(estimate = round(estimate, 4),
         across(c(p.value, p.fdr), \(x) signif(x, 3))) |>
  arrange(p.fdr)
if (nrow(sig) == 0) message("  (none)") else print(as.data.frame(sig))

# ── Forest plots ──────────────────────────────────────────────────────────────

save_forest_pair(
  all_res,
  out_dir = here("revisions", "figures", "epitope"),
  stem    = paste0("phyloglm_anychange_forest", sfx),
  title   = "Part 1: probability that any substitution occurred",
  xlab    = "Log-odds coefficient"
)

message("\nDone.")
