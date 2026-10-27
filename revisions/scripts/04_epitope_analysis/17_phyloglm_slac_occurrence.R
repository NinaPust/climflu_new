# revisions/scripts/04_epitope_analysis/17_phyloglm_slac_occurrence.R
#
# PURPOSE: Supplementary analysis to 16_phylolm_slac.R, retaining the branches
#          that 16 has to drop.
#
#          16 models log(dN) and log(dS), which are undefined when a branch
#          carries no change of that type — so it silently excludes 61% of
#          branches for the whole gene and 99% for the ridge. Reviewer 1
#          objected to exactly that exclusion. This script keeps every branch by
#          asking a question that a zero can answer:
#
#            any_N = did ANY non-synonymous change occur in this region?  (0/1)
#            any_S = did ANY synonymous change occur in this region?      (0/1)
#
#          Together the two scripts are a hurdle: 17 models whether change
#          happened, 16 models how fast it went given that it did.
#
# TIME IS REQUIRED HERE. Unlike omega — where N and S accumulate over the same
#   branch so duration cancels — a binary occurrence has nothing to cancel
#   against: a longer branch simply had more chances. Measured on H1,
#   P(N > 0) by branch-duration quartile:
#
#     full     27.1%  40.8%  46.3%  51.3%
#     HA1      16.3%  24.8%  27.8%  32.7%
#     HA2       7.9%  11.5%  15.2%  16.4%
#     epitope   2.9%   6.3%   8.2%   9.0%
#     ridge     0.7%   1.1%   1.5%   2.1%
#
#   Branches in the longest quartile are 2-3x more likely to show a change, in
#   every region. Since cor(scaled_temp, log t) = +0.11, omitting duration would
#   hand temperature a spurious positive coefficient in precisely the direction
#   the hypothesis predicts.
#
# EXPOSURE: log(EN * t) for any_N and log(ES * t) for any_S — the opportunity
#   for a non-synonymous change is the non-synonymous sites, and EN is itself
#   per branch (it drops where sequences are gappy). phyloglm has no offset()
#   support, so this enters as an ordinary estimated covariate rather than a
#   fixed offset. It is deliberately excluded from the FDR correction and from
#   the forest plots: it is a measurement construct, not a hypothesis. Its
#   fitted value is printed per model as a diagnostic.
#
# INTERPRETATION: a positive scaled_temp coefficient b means warmer locations
#   have exp(b) times the ODDS of carrying at least one such change, holding
#   exposure, the other climate variables and year fixed. One SD of temperature
#   is about 10 C in both genes. Example: b = 0.25 in the ridge any_N model
#   gives an odds ratio of 1.28, moving a 1.3% baseline occurrence to ~1.66%.
#   It is about OCCURRENCE, not amount — one change and three both read as 1.
#
# INPUT:
#   revisions/data/rates/{gene}_rates_slac.csv    from 15_slac_dnds.R
#   data/trees/timed/{gene}_timetree.nwk          branch durations (years)
#
# OUTPUT:
#   revisions/data/rates/phylolm/{gene}_phyloglm_occurrence_results.csv
#   revisions/data/rates/phylolm/all_phyloglm_occurrence_results.csv
#   revisions/figures/epitope/phyloglm_occurrence_forest_{N,S}.png  (+ _fdr)
#   revisions/figures/epitope/phyloglm_occurrence_forest_combined.png
#
# FLAGS:
#   --gene X        run a single gene (default: H1 H3)
#   --region X      run a single region (default: full HA1 HA2 epitope ridge)
#   --response X    any_N | any_S (default: both)
#   --no-year       drop as.factor(year). The ridge has 92-187 events against
#                   ~18 parameters, i.e. 5-7 per parameter against a rule of
#                   thumb of 10; dropping year leaves ~8 parameters and gets it
#                   to 11-15. Worth a second run for that region.
#   --skip-null     skip the year-only null fits. Halves the runtime but
#                   removes delta AIC, which is the main guard against reading
#                   too much into an individual p-value.
#   --min-events N  skip models with fewer than N events (default 20)
#   --threads N     BLAS/OMP cap (default 4; shared machine)
#   --alpha-bound X phyloglm log.alpha.bound (default 4, phyloglm's own). The
#                   first run pegged alpha at the resulting ceiling in every
#                   non-ridge model; try 6-8 to let it be estimated.
#   --plot-only     rebuild the figures from existing per-gene CSVs
#
# RUNTIME: this is the expensive script. phyloglm cost rises steeply with tip
#   count and every model here uses all tips (~10.9k H1, ~15.0k H3). Twenty
#   models, doubled to forty by the null fits. Run it in tmux. Per-gene results
#   are written after EACH model completes, so an interrupted run loses at most
#   one fit and can be resumed with --gene / --region / --response.

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

args       <- commandArgs(trailingOnly = TRUE)
genes      <- if ("--gene"     %in% args) args[which(args == "--gene")     + 1] else c("H1", "H3")
regions    <- if ("--region"   %in% args) args[which(args == "--region")   + 1] else
              c("full", "HA1", "HA2", "epitope", "ridge")
responses  <- if ("--response" %in% args) args[which(args == "--response") + 1] else
              c("any_N", "any_S")
use_year   <- !("--no-year"   %in% args)
fit_null   <- !("--skip-null" %in% args)
plot_only  <- "--plot-only"   %in% args
min_events <- if ("--min-events" %in% args) as.integer(args[which(args == "--min-events") + 1]) else 20L
n_thread   <- if ("--threads"    %in% args) as.integer(args[which(args == "--threads")    + 1]) else 4L
# phyloglm's default of 4 pegged alpha in every non-ridge model of the first
# run; raising it lets alpha be estimated rather than pinned at the ceiling.
alpha_bound <- if ("--alpha-bound" %in% args) as.numeric(args[which(args == "--alpha-bound") + 1]) else 4

cap_blas_threads(n_thread)

PRED <- c("scaled_temp", "scaled_humidity",
          "scaled_delta_temp", "scaled_delta_humidity",
          "scaled_sampling_density", "scaled_pop_density")

RESPONSE_LABELS <- c(
  any_N = "Occurrence of a non-synonymous change",
  any_S = "Occurrence of a synonymous change"
)

# Which SLAC columns each response is built from: count, and the matching
# expected-site count that forms its exposure.
RESPONSE_SPEC <- list(
  any_N = list(count = "N", sites = "EN"),
  any_S = list(count = "S", sites = "ES")
)

out_dir <- PATHS$rates_phylolm
fig_dir <- here("revisions", "figures", "epitope")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_file <- function(g)
  file.path(out_dir, paste0(g, "_phyloglm_occurrence_results.csv"))

# ── Fit ───────────────────────────────────────────────────────────────────────

if (!plot_only) for (gene in genes) {
  message("\n========== ", gene, " ==========")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_slac.csv"))
  tree_path  <- file.path(here(), "data", "trees", "timed",
                          paste0(gene, "_timetree.nwk"))
  if (!file.exists(rates_path)) {
    warning(gene, ": ", basename(rates_path), " not found — run 15_slac_dnds.R")
    next
  }
  if (!file.exists(tree_path)) { warning(gene, ": timetree not found"); next }

  df        <- suppressWarnings(read_csv(rates_path, show_col_types = FALSE,
                                         guess_max = 50000))
  full_tree <- read.tree(tree_path)
  # Durations from the UNPRUNED tree: pruning reassigns terminal branch lengths,
  # so exposure must be measured before any model-specific subsetting.
  df$t <- terminal_branch_times(full_tree)[df$accession]

  gene_rows <- list()

  for (r in regions) {
    for (resp in responses) {
      spec     <- RESPONSE_SPEC[[resp]]
      cnt_col  <- paste0("slac_", r, "_", spec$count)
      site_col <- paste0("slac_", r, "_", spec$sites)
      if (!all(c(cnt_col, site_col) %in% names(df))) {
        message("  [", r, " ", resp, "] SLAC columns missing — skipping"); next
      }

      sub <- df
      sub$response     <- as.integer(sub[[cnt_col]] > 0)
      sub$log_exposure <- log(sub[[site_col]] * sub$t)

      keep <- intersect(c("accession", "response", "log_exposure", PRED,
                          if (use_year) "year"), names(sub))
      sub  <- sub[, keep, drop = FALSE]
      # EN/ES can be 0 where a branch has no scoreable codon in the region, and
      # t can be 0; both make log_exposure non-finite.
      sub  <- sub[stats::complete.cases(sub) & is.finite(sub$log_exposure), ,
                  drop = FALSE]
      if (nrow(sub) < 100) {
        message("  [", r, " ", resp, "] <100 usable branches — skipping"); next
      }

      prepped  <- prep_tree_data(sub, full_tree, quiet = TRUE)
      sub_df   <- prepped$data
      sub_tree <- prepped$tree

      n_ev <- sum(sub_df$response == 1)
      if (n_ev < min_events || n_ev > nrow(sub_df) - min_events) {
        message(sprintf("  [%-7s %-5s] n = %5d  events = %4d  — too few (or too many) to fit, skipping",
                        r, resp, nrow(sub_df), n_ev))
        next
      }

      preds     <- intersect(PRED, names(sub_df))
      year_term <- if (use_year && "year" %in% names(sub_df)) "as.factor(year)" else NULL
      n_par     <- length(preds) + 2L +
                   (if (!is.null(year_term)) length(unique(sub_df$year)) - 1L else 0L)
      fml <- as.formula(paste("response ~",
                              paste(c(preds, "log_exposure", year_term),
                                    collapse = " + ")))

      message(sprintf("  [%-7s %-5s] n = %5d  events = %4d (%.1f%%)  ev/par = %4.1f  fitting...",
                      r, resp, nrow(sub_df), n_ev, 100 * n_ev / nrow(sub_df),
                      n_ev / n_par))
      if (n_ev / n_par < 10)
        message("      NOTE: under 10 events per parameter — estimates ",
                "unstable here; consider --no-year for this region")

      # The null keeps log_exposure and the year factor and drops ONLY the
      # climate block, so delta AIC answers "what do the climate predictors
      # add". A year-only null would also be penalised for losing the exposure
      # term — which strongly predicts occurrence — and would inflate delta AIC
      # for reasons having nothing to do with climate.
      t0  <- proc.time()
      res <- tryCatch(
        run_phyloglm(formula = fml, data = sub_df, tree = sub_tree,
                     subset = r, response = resp, gene = gene,
                     fit_null = fit_null, null_terms = "log_exposure",
                     log_alpha_bound = alpha_bound),
        error = function(e) { message("      FAILED: ", conditionMessage(e)); NULL }
      )
      if (is.null(res)) next
      el <- round((proc.time() - t0)[["elapsed"]] / 60, 1)

      tmp  <- res[res$term == "scaled_temp", ]
      expc <- res$estimate[res$term == "log_exposure"]
      dAIC <- if (!is.na(res$AIC_null[1])) res$AIC_null[1] - res$AIC[1] else NA_real_
      message(sprintf("      temp = %+.4f (p = %-9.3g)  exposure = %+.3f  dAIC = %s  [%s min]",
                      tmp$estimate, tmp$p.value, expc,
                      if (is.na(dAIC)) "NA" else sprintf("%.1f", dAIC), el))

      res$obs_per_par <- nrow(sub_df) / n_par
      gene_rows[[paste(r, resp)]] <- res

      # Written after every model: this run is long, and an interruption should
      # cost at most one fit.
      write_csv(bind_rows(gene_rows), gene_file(gene))
    }
  }

  if (length(gene_rows) == 0) {
    message("  no successful fits for ", gene)
  } else {
    message("  Written: ", basename(gene_file(gene)), "  (",
            length(gene_rows), " models)")
  }
}

# ── Combine and FDR-correct ───────────────────────────────────────────────────

all_res <- lapply(genes, function(g) {
  p <- gene_file(g)
  if (file.exists(p)) suppressWarnings(read_csv(p, show_col_types = FALSE)) else NULL
}) |> bind_rows()

if (nrow(all_res) == 0) stop("No results found. Run without --plot-only first.")

all_res$p.fdr <- NA_real_
pr <- all_res$term %in% PRED
all_res$p.fdr[pr] <- p.adjust(all_res$p.value[pr], method = "BH")

combined <- file.path(out_dir, "all_phyloglm_occurrence_results.csv")
write_csv(all_res, combined)
message("\nWritten: ", basename(combined))

# ── Report ────────────────────────────────────────────────────────────────────

message("\nTemperature coefficient, all models (OR = odds ratio per 1 SD, ~10 C):")
all_res |>
  filter(term == "scaled_temp") |>
  mutate(estimate = round(estimate, 4), OR = round(exp(estimate), 3),
         p = signif(p.value, 3), fdr = signif(p.fdr, 3),
         dAIC = ifelse(is.na(AIC_null), NA, round(AIC_null - AIC, 1)),
         sig = ifelse(!is.na(p.fdr) & p.fdr < 0.05, "*", "")) |>
  select(gene, region = subset, response, n, n_events, estimate, OR,
         p, fdr, sig, dAIC) |>
  arrange(response, gene, region) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nExposure coefficients (opportunity term; not a hypothesis):")
all_res |>
  filter(term == "log_exposure") |>
  mutate(estimate = round(estimate, 3), se = round(std_error, 3)) |>
  select(gene, region = subset, response, estimate, se, n) |>
  arrange(response, gene, region) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nClimate terms significant after BH-FDR:")
sig <- all_res |>
  filter(term %in% PRED, !is.na(p.fdr), p.fdr < 0.05) |>
  mutate(estimate = round(estimate, 4), OR = round(exp(estimate), 3),
         p.fdr = signif(p.fdr, 3)) |>
  select(gene, region = subset, response, term, estimate, OR, p.fdr,
         n, n_events) |>
  arrange(p.fdr)
if (nrow(sig) == 0) message("  (none)") else print(as.data.frame(sig), row.names = FALSE)

# ── Forest plots ──────────────────────────────────────────────────────────────

#' Prepare results for plotting: CIs, labels, significance flags.
prep_plot <- function(res) {
  res |>
    filter(term %in% PRED) |>
    mutate(
      ci_lo     = estimate - 1.96 * std_error,
      ci_hi     = estimate + 1.96 * std_error,
      sig_raw   = !is.na(p.value) & p.value < 0.05,
      sig_fdr   = !is.na(p.fdr)   & p.fdr   < 0.05,
      predictor = factor(PREDICTOR_LABELS[term], levels = PREDICTOR_LABELS),
      region    = factor(REGION_LABELS[subset], levels = REGION_LEVELS),
      subtype   = ifelse(gene %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2")
    ) |>
    filter(!is.na(region))
}

#' Forest plot for one response, faceted by gene, with delta-AIC annotations.
forest_one <- function(res, resp, sig_col) {
  d <- prep_plot(res) |> filter(response == resp)
  if (!nrow(d)) return(NULL)
  d$point_fill <- ifelse(d[[sig_col]],
                         PREDICTOR_COLOURS[as.character(d$predictor)], "white")
  lab <- d |>
    distinct(subtype, region, AIC, AIC_null) |>
    filter(!is.na(AIC_null)) |>
    mutate(txt = sprintf("ΔAIC = %.1f", AIC_null - AIC))

  p <- ggplot(d, aes(x = estimate, y = region,
                     colour = predictor, group = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.7,
                   position = position_dodge(width = 0.7)) +
    geom_point(aes(fill = point_fill), shape = 21, size = 2.6, stroke = 0.8,
               position = position_dodge(width = 0.7))
  if (nrow(lab) > 0)
    p <- p + geom_text(data = lab, aes(x = Inf, y = region, label = txt),
                       hjust = 1.05, inherit.aes = FALSE, size = 2.5,
                       colour = "grey40")
  p +
    scale_colour_manual(values = PREDICTOR_COLOURS, name = NULL) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev, drop = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.34))) +
    facet_wrap(~ subtype) +
    labs(title = RESPONSE_LABELS[[resp]],
         subtitle = if (sig_col == "sig_fdr")
           "Log-odds coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05"
         else "Log-odds coefficients ± 95% CI  |  Filled = p < 0.05",
         x = "Log-odds coefficient", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold"),
          strip.background = element_blank(),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    guides(colour = guide_legend(nrow = 2))
}

for (resp in intersect(responses, unique(all_res$response))) {
  for (v in list(list(col = "sig_raw", sfx = ""),
                 list(col = "sig_fdr", sfx = "_fdr"))) {
    p <- forest_one(all_res, resp, v$col)
    if (is.null(p)) next
    stem <- paste0("phyloglm_occurrence_forest_", sub("^any_", "", resp), v$sfx)
    ggsave(file.path(fig_dir, paste0(stem, ".png")), p,
           width = 9, height = 5, dpi = 200)
    message("  Saved: ", stem, ".png")
  }
}

# Combined panel: response (rows) x gene (columns)
comb <- prep_plot(all_res) |>
  mutate(resp_lab = factor(RESPONSE_LABELS[response], levels = RESPONSE_LABELS))
comb$point_fill <- ifelse(comb$sig_fdr,
                          PREDICTOR_COLOURS[as.character(comb$predictor)], "white")
comb_lab <- comb |>
  distinct(resp_lab, subtype, region, AIC, AIC_null) |>
  filter(!is.na(AIC_null)) |>
  mutate(txt = sprintf("ΔAIC = %.1f", AIC_null - AIC))

p_comb <- ggplot(comb, aes(x = estimate, y = region,
                           colour = predictor, group = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.6,
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = point_fill), shape = 21, size = 2.2, stroke = 0.7,
             position = position_dodge(width = 0.7))
if (nrow(comb_lab) > 0)
  p_comb <- p_comb +
    geom_text(data = comb_lab, aes(x = Inf, y = region, label = txt),
              hjust = 1.05, inherit.aes = FALSE, size = 2.1, colour = "grey40")
p_comb <- p_comb +
  scale_colour_manual(values = PREDICTOR_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.34))) +
  facet_grid(resp_lab ~ subtype, scales = "free_x") +
  labs(title = "Occurrence of substitution vs climate (zeroes retained)",
       subtitle = "Log-odds coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05",
       x = "Log-odds coefficient", y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  guides(colour = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "phyloglm_occurrence_forest_combined.png"), p_comb,
       width = 9, height = 7, dpi = 200)
message("  Saved: phyloglm_occurrence_forest_combined.png")

message("\nDone.")
