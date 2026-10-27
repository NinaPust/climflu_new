# revisions/scripts/04_epitope_analysis/12b_phyloglm_ns_given_change.R
#
# HURDLE PART 2 — "given that a substitution occurred, was it protein-changing?"
#
# PURPOSE: Phylogenetic logistic regression of substitution TYPE against climate
#          predictors. This is the SELECTION half of the hurdle, and it is a
#          model of omega.
#
#          Population: terminal branches carrying EXACTLY ONE substitution.
#          Response:   1 if that substitution was non-synonymous, 0 if synonymous.
#
# WHY THIS IS TIME-FREE — the point of the whole exercise:
#          Both outcomes required the same branch to exist. The branch already
#          "spent" its duration acquiring one substitution; whether that
#          substitution was NS or S does not depend on how long the branch was.
#          So no exposure term is needed and none is included. This is what
#          escapes the artefact in which log(count / branch_time) collapses to
#          -log(branch_time) plus noise, because counts on terminal branches are
#          almost always 0 or 1 and the Nei-Gojobori denominator is near-constant
#          within a region.
#
# WHY EXACTLY ONE SUBSTITUTION:
#          With one substitution, "was it NS?" is the COMPLETE data for that
#          branch — the Bernoulli is exact, not an approximation of a binomial.
#          Collapsing a branch with 1 NS and 3 S to "yes" would bias omega
#          upward. Restricting costs little in these regions (~79-99% of
#          informative tips carry exactly one substitution) but would cost ~45%
#          for the whole gene, which is why `full` is excluded here and left to
#          Part 1 and to the count model in 12c.
#
# WHY NOT A CONTINUOUS PER-TIP OMEGA:
#          omega = (Nd/N)/(Sd/S) needs at least one substitution of EACH type on
#          the same branch. Only 0.4-2% of tips qualify in these regions (H3
#          ridge: zero tips). Nei-Gojobori normalisation does not rescue this —
#          the zero sits in the numerator of pS = Sd/S, and dividing zero by a
#          nonzero site count is still zero. A tip with one substitution carries
#          one bit of information; a continuous ratio cannot be extracted from
#          one bit. This model is honest about that: one bit per tip, pooled
#          across tips, omega estimated at the population level.
#
# READING THE INTERCEPT AS OMEGA:
#          log-odds(NS) - log(N/S) = log(omega), where N/S is the Nei-Gojobori
#          expected site ratio. Because N/S is near-constant within a region it
#          cannot enter as a per-tip offset (and phyloglm has no offset()
#          support anyway) — it simply shifts the intercept. So:
#              baseline log(omega) = intercept - mean(log(N/S))
#          is reported per gene x region below and stored in the results CSV.
#          Crucially, a constant shift does NOT affect the slopes: the climate
#          coefficients read directly as effects on log(omega).
#
# INTERPRETATION: a positive temperature coefficient means omega is HIGHER at
#          warmer temperatures — purifying selection is more relaxed, or
#          diversifying selection stronger. The sign is meaningful. Exponentiate
#          for effect size: exp(b) is the multiplicative change in omega per one
#          standard deviation of the predictor (1 SD ~ 10 C for both genes).
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/phylolm/{gene}_phyloglm_nsgiven_results.csv
#   data/rates/phylolm/all_phyloglm_nsgiven_results.csv
#   figures/epitope/phyloglm_nsgiven_forest.png
#   figures/epitope/phyloglm_nsgiven_forest_fdr.png
#
# FLAGS:
#   --gene X      run a single gene (default: H1 H3)
#   --region X    run a single region (default: HA2 epitope ridge)
#   --allow-full  also run the whole gene, accepting the exactly-one
#                 restriction's ~45% data loss there (off by default)
#   --any-ns      relax the population to ALL tips with >=1 substitution, with
#                 response "was any of it NS". Uses more tips but biases omega
#                 upward on multi-substitution branches. Diagnostic only.
#   --no-year     drop as.factor(year); advisable for the ridge (~200-260 obs)
#   --fit-null    also fit a year-only model for delta-AIC (doubles runtime)
#   --threads N   BLAS/OMP cap (default 4; shared machine)
#
# RUNTIME: much cheaper than 12a — n is 1.2k-3.6k rather than every tip.

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
            c("HA2", "epitope", "ridge")
if ("--allow-full" %in% args && !("--region" %in% args))
  regions <- c("full", regions)
any_ns   <- "--any-ns"    %in% args
use_year <- !("--no-year" %in% args)
fit_null <- "--fit-null"  %in% args
n_thread <- if ("--threads" %in% args) as.integer(args[which(args == "--threads") + 1]) else 4L

cap_blas_threads(n_thread)

if (any_ns)
  message("MODE: --any-ns — population is tips with >=1 substitution, ",
          "response is 'any NS'. Biased upward on multi-substitution branches; ",
          "diagnostic only.")

predictor_terms <- c(
  "scaled_temp", "scaled_humidity",
  "scaled_delta_temp", "scaled_delta_humidity",
  "scaled_sampling_density", "scaled_pop_density"
)

out_dir <- file.path(PATHS$rates_phylolm)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Fit per gene x region ─────────────────────────────────────────────────────

for (gene in genes) {
  message("\n[", gene, "] phyloglm — NS given change (hurdle part 2)")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
  if (!file.exists(rates_path)) { warning("Rates not found: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Tree not found: ",  tree_path);  next }

  df        <- read_csv(rates_path, show_col_types = FALSE)
  full_tree <- read.tree(tree_path)

  gene_results <- list()

  for (r in regions) {
    ns_col <- paste0(r, "_ns_count"); s_col <- paste0(r, "_s_count")
    N_col  <- paste0(r, "_N_sites");  S_col <- paste0(r, "_S_sites")
    needed <- c(ns_col, s_col, N_col, S_col)
    if (!all(needed %in% names(df))) {
      message("  [", r, "] missing columns (",
              paste(setdiff(needed, names(df)), collapse = ", "),
              ") — re-run 09_substitutions.R; skipping")
      next
    }

    sub <- df
    sub$tot_sub <- sub[[ns_col]] + sub[[s_col]]

    # Population and response depend on mode.
    if (any_ns) {
      sub <- sub[!is.na(sub$tot_sub) & sub$tot_sub >= 1, , drop = FALSE]
      sub$ns_given <- as.integer(sub[[ns_col]] > 0)
    } else {
      sub <- sub[!is.na(sub$tot_sub) & sub$tot_sub == 1, , drop = FALSE]
      sub$ns_given <- as.integer(sub[[ns_col]] == 1)
    }

    # Nei-Gojobori site ratio, for converting the intercept to baseline omega.
    sub$log_NS_ratio <- log(sub[[N_col]] / sub[[S_col]])

    keep <- c("accession", "ns_given", "log_NS_ratio", predictor_terms,
              if (use_year) "year")
    keep <- intersect(keep, names(sub))
    sub  <- sub[, keep, drop = FALSE]
    sub  <- sub[stats::complete.cases(sub) &
                is.finite(sub$log_NS_ratio), , drop = FALSE]

    if (nrow(sub) < 50) { message("  [", r, "] <50 usable tips, skipping"); next }

    prepped  <- prep_tree_data(sub, full_tree)
    sub_df   <- prepped$data
    sub_tree <- prepped$tree

    n_ns  <- sum(sub_df$ns_given == 1)
    pct   <- 100 * n_ns / nrow(sub_df)
    mlr   <- mean(sub_df$log_NS_ratio)
    message("  [", r, "] n = ", nrow(sub_df), "  NS = ", n_ns,
            " (", round(pct, 1), "%)  mean N/S = ", round(exp(mlr), 2))

    if (n_ns < 20 || n_ns > nrow(sub_df) - 20) {
      message("    too few NS (or S) outcomes for a stable fit, skipping")
      next
    }

    preds     <- intersect(predictor_terms, names(sub_df))
    year_term <- if (use_year && "year" %in% names(sub_df)) "as.factor(year)" else NULL
    n_par     <- length(preds) + 1L + (if (!is.null(year_term))
                                         length(unique(sub_df$year)) - 1L else 0L)
    n_min     <- min(n_ns, nrow(sub_df) - n_ns)
    if (n_min / n_par < 10)
      message("    NOTE: ", round(n_min / n_par, 1), " minority outcomes per ",
              "parameter (<10) — unstable; consider --no-year")

    # No exposure term: see header. Branch time cancels between the two outcomes.
    fml <- as.formula(paste("ns_given ~",
                            paste(c(preds, year_term), collapse = " + ")))

    t0  <- proc.time()
    res <- tryCatch(
      run_phyloglm(formula = fml, data = sub_df, tree = sub_tree,
                   subset = r, response = "ns_given", gene = gene,
                   fit_null = fit_null),
      error = function(e) { message("    FAILED: ", conditionMessage(e)); NULL }
    )
    if (is.null(res)) next
    message("    fitted in ", round((proc.time() - t0)[["elapsed"]]), "s")

    # Baseline omega from the intercept, shifted by the (near-constant) site ratio.
    icpt <- res$estimate[res$term == "(Intercept)"]
    res$mean_log_NS_ratio  <- mlr
    res$baseline_log_omega <- if (length(icpt) == 1) icpt - mlr else NA_real_
    res$baseline_omega     <- exp(res$baseline_log_omega)
    if (length(icpt) == 1)
      message("    baseline log(omega) = ", signif(icpt, 4), " - ",
              signif(mlr, 4), " = ", signif(icpt - mlr, 4),
              "  ->  omega = ", signif(exp(icpt - mlr), 4))

    gene_results[[r]] <- res
  }

  if (length(gene_results) == 0) { message("  no successful fits"); next }
  gene_df  <- bind_rows(gene_results)
  gene_out <- file.path(out_dir, paste0(gene, "_phyloglm_nsgiven_results.csv"))
  write_csv(gene_df, gene_out)
  message("  Written: ", basename(gene_out))
}

# ── Combine, FDR-correct, report ──────────────────────────────────────────────

all_res <- lapply(genes, function(g) {
  p <- file.path(out_dir, paste0(g, "_phyloglm_nsgiven_results.csv"))
  if (file.exists(p)) read_csv(p, show_col_types = FALSE) else NULL
}) |> bind_rows()

if (nrow(all_res) == 0) stop("No results produced.")

all_res <- add_fdr(all_res, predictor_terms)
combined <- file.path(out_dir, "all_phyloglm_nsgiven_results.csv")
write_csv(all_res, combined)
message("\nWritten: ", basename(combined))

message("\nBaseline omega per gene x region:")
all_res |>
  distinct(gene, subset, baseline_omega, n, n_events) |>
  mutate(baseline_omega = signif(baseline_omega, 3)) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nSignificant after BH-FDR (p.fdr < 0.05):")
sig <- all_res |>
  filter(term %in% predictor_terms, !is.na(p.fdr), p.fdr < 0.05) |>
  mutate(omega_multiplier = round(exp(estimate), 3),
         estimate = round(estimate, 4),
         across(c(p.value, p.fdr), \(x) signif(x, 3))) |>
  select(gene, subset, term, estimate, omega_multiplier, p.value, p.fdr, n) |>
  arrange(p.fdr)
if (nrow(sig) == 0) message("  (none)") else print(as.data.frame(sig))

# ── Forest plots ──────────────────────────────────────────────────────────────

save_forest_pair(
  all_res,
  out_dir = here("revisions", "figures", "epitope"),
  stem    = "phyloglm_nsgiven_forest",
  title   = "Part 2: probability a substitution was non-synonymous (omega)",
  xlab    = "Coefficient on log(omega)"
)

message("\nDone.")
