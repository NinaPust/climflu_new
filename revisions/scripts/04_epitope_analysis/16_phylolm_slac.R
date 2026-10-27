# revisions/scripts/04_epitope_analysis/16_phylolm_slac.R
#
# PURPOSE: Phylogenetic regression of SLAC-derived dN, dS and omega against
#          climate, per gene and per HA region. This is the primary analysis
#          behind the revised figure 1D.
#
#          Inputs come from 15_slac_dnds.R, so the substitution counts are
#          maximum-likelihood (HyPhy SLAC) rather than parsimony, and the
#          expected site counts (EN, ES) are per branch, from the fitted codon
#          model, restricted to codons scoreable on that branch.
#
# THE THREE RESPONSES:
#
#   log_dN_rate = log( N / (EN * t) )      non-syn substitutions per non-syn
#                                          site per year     [requires N > 0]
#   log_dS_rate = log( S / (ES * t) )      synonymous equivalent  [S > 0]
#   log_omega   = log( (N/EN) / (S/ES) )   selection intensity  [N > 0 AND S > 0]
#
#   t is the terminal branch duration in years from the timetree. SLAC itself
#   has no notion of calendar time — its branch lengths are in substitutions per
#   site — so time has to be brought in here. Verified on H1: mean N rises from
#   0.32 to 0.83 across branch-duration quartiles, so duration genuinely has to
#   be accounted for in a rate.
#
#   omega needs no time term: N and S accumulated over the same branch, so t
#   cancels algebraically. That is a property of the quantity, not a choice.
#
# ZEROES ARE DROPPED. A branch with no non-synonymous change has dN = 0 and
#   log(dN) undefined, so it leaves the dN model; likewise dS, and omega needs
#   both. This is deliberate and matches the submitted analysis: the question
#   here is what drives the rate among lineages that did change. A supplementary
#   analysis retaining the zeroes is handled separately. Reported n per model
#   makes the size of the exclusion explicit.
#
# KNOWN PROPERTY, worth reading results against: among branches with N > 0, the
#   fraction carrying exactly N = 1 runs from 76% (whole gene) to 95-99%
#   (epitope, ridge). For those, log(N) = 0 and the response reduces to
#   -log(EN) - log(t). EN varies genuinely per branch (missing data), but branch
#   duration still drives much of the variance in the short regions. Delta AIC
#   is the check: a near-zero value means the climate terms add nothing over the
#   year factor, whatever the individual p-values look like.
#
# MODEL (matching scripts/03_rates_analysis/02_phylolm.R):
#   full:  response ~ 6 scaled climate predictors + as.factor(year)
#   null:  response ~ as.factor(year)            -> delta AIC = AIC_null - AIC
#   phylolm with Pagel's lambda; tree pruned per model to branches with a
#   defined response; zero-length terminal branches dropped.
#
# FDR: Benjamini-Hochberg across every climate term, over all genes x regions
#      x responses.
#
# INPUT:
#   revisions/data/rates/{gene}_rates_slac.csv    from 15_slac_dnds.R
#   data/trees/timed/{gene}_timetree.nwk          branch durations (years)
#
# OUTPUT:
#   revisions/data/rates/phylolm/{gene}_phylolm_slac_results.csv
#   revisions/data/rates/phylolm/all_phylolm_slac_results.csv
#   revisions/figures/epitope/phylolm_slac_forest_{dN,dS,omega}.png  (+ _fdr)
#   revisions/figures/epitope/phylolm_slac_forest_combined.png
#
# FLAGS:
#   --gene X        run a single gene (default: H1 H3)
#   --region X      run a single region (default: full HA1 HA2 epitope ridge)
#   --response X    log_dN_rate | log_dS_rate | log_omega (default: all three)
#   --no-year       drop as.factor(year) (frees ~10 parameters; worth trying
#                   for the ridge, where n is 92-187)
#   --min-n N       skip models with fewer than N observations (default 50)
#   --threads N     BLAS/OMP cap (default 4; shared machine)
#
# RUNTIME: 30 models plus 30 null fits; a few minutes at these sample sizes.

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(viridisLite)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "phylolm_utils.R"))
source(here("revisions", "R", "hurdle_utils.R"))
source(here("revisions", "R", "forest_utils.R"))

args     <- commandArgs(trailingOnly = TRUE)
genes    <- if ("--gene"     %in% args) args[which(args == "--gene")     + 1] else c("H1", "H3")
regions  <- if ("--region"   %in% args) args[which(args == "--region")   + 1] else
            c("full", "HA1", "HA2", "epitope", "ridge")
responses <- if ("--response" %in% args) args[which(args == "--response") + 1] else
            c("log_dN_rate", "log_dS_rate", "log_omega")
use_year <- !("--no-year" %in% args)
min_n    <- if ("--min-n"   %in% args) as.integer(args[which(args == "--min-n")   + 1]) else 50L
n_thread <- if ("--threads" %in% args) as.integer(args[which(args == "--threads") + 1]) else 4L

cap_blas_threads(n_thread)

PRED <- c("scaled_temp", "scaled_humidity",
          "scaled_delta_temp", "scaled_delta_humidity",
          "scaled_sampling_density", "scaled_pop_density")

RESPONSE_LABELS <- c(
  log_dN_rate = "dN  (non-synonymous rate)",
  log_dS_rate = "dS  (synonymous rate)",
  log_omega   = "omega  (dN/dS)"
)

out_dir <- PATHS$rates_phylolm
fig_dir <- here("revisions", "figures", "epitope")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


#' Build the three responses for one region.
#'
#' @param d  Data frame with slac_{region}_{N,S,EN,ES} columns and `t` (years).
#' @param r  Region key.
#' @return `d` with log_dN_rate, log_dS_rate and log_omega added.
add_responses <- function(d, r) {
  N  <- d[[paste0("slac_", r, "_N")]];  S  <- d[[paste0("slac_", r, "_S")]]
  EN <- d[[paste0("slac_", r, "_EN")]]; ES <- d[[paste0("slac_", r, "_ES")]]

  okN <- !is.na(N) & !is.na(EN) & EN > 0 & !is.na(d$t) & d$t > 0 & N > 0
  okS <- !is.na(S) & !is.na(ES) & ES > 0 & !is.na(d$t) & d$t > 0 & S > 0
  okO <- !is.na(N) & !is.na(S) & !is.na(EN) & !is.na(ES) &
         EN > 0 & ES > 0 & N > 0 & S > 0

  d$log_dN_rate <- ifelse(okN, log(N / (EN * d$t)), NA_real_)
  d$log_dS_rate <- ifelse(okS, log(S / (ES * d$t)), NA_real_)
  # No t: it cancels between numerator and denominator.
  d$log_omega   <- ifelse(okO, log((N / EN) / (S / ES)), NA_real_)
  d
}

# ── Fit ───────────────────────────────────────────────────────────────────────

all_results <- list()

for (gene in genes) {
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
  df$t      <- terminal_branch_times(full_tree)[df$accession]

  gene_rows <- list()

  for (r in regions) {
    need <- paste0("slac_", r, c("_N", "_S", "_EN", "_ES"))
    if (!all(need %in% names(df))) {
      message("  [", r, "] SLAC columns missing — skipping"); next
    }
    dr <- add_responses(df, r)

    for (resp in responses) {
      keep <- intersect(c("accession", resp, PRED, if (use_year) "year"),
                        names(dr))
      sub  <- dr[, keep, drop = FALSE]
      sub  <- sub[stats::complete.cases(sub) & is.finite(sub[[resp]]), ,
                  drop = FALSE]

      if (nrow(sub) < min_n) {
        message(sprintf("  [%-7s %-11s] n = %d  (< %d, skipped)",
                        r, resp, nrow(sub), min_n))
        next
      }

      prepped  <- prep_tree_data(sub, full_tree, quiet = TRUE)
      sub_df   <- prepped$data
      sub_tree <- prepped$tree

      preds     <- intersect(PRED, names(sub_df))
      year_term <- if (use_year && "year" %in% names(sub_df)) "as.factor(year)" else NULL
      n_par     <- length(preds) + 1L +
                   (if (!is.null(year_term)) length(unique(sub_df$year)) - 1L else 0L)
      fml <- as.formula(paste(resp, "~",
                              paste(c(preds, year_term), collapse = " + ")))

      t0  <- proc.time()
      res <- tryCatch(
        run_phylolm(formula = fml, data = sub_df, tree = sub_tree,
                    method = PARAMS$phylolm_method,
                    subset = r, response = resp, gene = gene),
        error = function(e) { message("    FAILED: ", conditionMessage(e)); NULL }
      )
      if (is.null(res)) next
      el <- round((proc.time() - t0)[["elapsed"]])

      tmp  <- res[res$term == "scaled_temp", ]
      dAIC <- res$AIC_null[1] - res$AIC[1]
      message(sprintf("  [%-7s %-11s] n = %5d  temp = %+.4f (p = %-9.3g) dAIC = %7.1f  obs/par = %4.1f  [%ds]",
                      r, resp, nrow(sub_df), tmp$estimate, tmp$p.value,
                      dAIC, nrow(sub_df) / n_par, el))
      if (nrow(sub_df) / n_par < 10)
        message("      NOTE: fewer than 10 observations per parameter — ",
                "estimates unstable; consider --no-year for this region")

      res$obs_per_par <- nrow(sub_df) / n_par
      gene_rows[[paste(r, resp)]] <- res
    }
  }

  if (length(gene_rows) == 0) { message("  no successful fits"); next }
  gdf <- bind_rows(gene_rows)
  gout <- file.path(out_dir, paste0(gene, "_phylolm_slac_results.csv"))
  write_csv(gdf, gout)
  message("  Written: ", basename(gout))
  all_results[[gene]] <- gdf
}

if (length(all_results) == 0) stop("No models fitted.")

all_res <- bind_rows(all_results)
all_res$p.fdr <- NA_real_
pr <- all_res$term %in% PRED
all_res$p.fdr[pr] <- p.adjust(all_res$p.value[pr], method = "BH")

combined <- file.path(out_dir, "all_phylolm_slac_results.csv")
write_csv(all_res, combined)
message("\nWritten: ", basename(combined))

# ── Report ────────────────────────────────────────────────────────────────────

message("\nTemperature coefficient, all models:")
all_res |>
  filter(term == "scaled_temp") |>
  mutate(estimate = round(estimate, 4), p = signif(p.value, 3),
         fdr = signif(p.fdr, 3), dAIC = round(AIC_null - AIC, 1),
         sig = ifelse(!is.na(p.fdr) & p.fdr < 0.05, "*", "")) |>
  select(gene, region = subset, response, n, estimate, p, fdr, sig, dAIC) |>
  arrange(response, gene, region) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nAll climate terms significant after BH-FDR:")
sig <- all_res |>
  filter(term %in% PRED, !is.na(p.fdr), p.fdr < 0.05) |>
  mutate(estimate = round(estimate, 4), p.fdr = signif(p.fdr, 3)) |>
  select(gene, region = subset, response, term, estimate, p.fdr, n) |>
  arrange(p.fdr)
if (nrow(sig) == 0) message("  (none)") else print(as.data.frame(sig), row.names = FALSE)

# ── Forest plots ──────────────────────────────────────────────────────────────

#' Forest plot for one response, faceted by gene, with delta-AIC annotations.
forest_one <- function(res, resp, sig_col) {
  d <- res |> filter(response == resp, term %in% PRED)
  if (!nrow(d)) return(NULL)
  d <- d |>
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
  d$point_fill <- ifelse(d[[sig_col]],
                         PREDICTOR_COLOURS[as.character(d$predictor)], "white")

  # Placed at the panel's right edge rather than a computed x, so it sits
  # correctly whatever the coefficient scale; the x expansion makes room.
  lab <- d |>
    distinct(subtype, region, AIC, AIC_null) |>
    mutate(txt = sprintf("ΔAIC = %.1f", AIC_null - AIC))

  ggplot(d, aes(x = estimate, y = region, colour = predictor, group = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.7,
                   position = position_dodge(width = 0.7)) +
    geom_point(aes(fill = point_fill), shape = 21, size = 2.6, stroke = 0.8,
               position = position_dodge(width = 0.7)) +
    geom_text(data = lab, aes(x = Inf, y = region, label = txt),
              hjust = 1.05, inherit.aes = FALSE, size = 2.5, colour = "grey40") +
    scale_colour_manual(values = PREDICTOR_COLOURS, name = NULL) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev, drop = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.05, 0.28))) +
    facet_wrap(~ subtype) +
    labs(title = RESPONSE_LABELS[[resp]],
         subtitle = if (sig_col == "sig_fdr")
           "Coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05"
         else "Coefficients ± 95% CI  |  Filled = p < 0.05",
         x = "Regression coefficient", y = NULL) +
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
    stem <- paste0("phylolm_slac_forest_",
                   sub("^log_", "", sub("_rate$", "", resp)), v$sfx)
    ggsave(file.path(fig_dir, paste0(stem, ".png")), p,
           width = 9, height = 5, dpi = 200)
    message("  Saved: ", stem, ".png")
  }
}

# Combined panel: response (rows) x gene (columns)
comb <- all_res |>
  filter(term %in% PRED) |>
  mutate(ci_lo = estimate - 1.96 * std_error,
         ci_hi = estimate + 1.96 * std_error,
         sig_fdr = !is.na(p.fdr) & p.fdr < 0.05,
         predictor = factor(PREDICTOR_LABELS[term], levels = PREDICTOR_LABELS),
         region    = factor(REGION_LABELS[subset], levels = REGION_LEVELS),
         subtype   = ifelse(gene %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2"),
         resp_lab  = factor(RESPONSE_LABELS[response],
                            levels = RESPONSE_LABELS)) |>
  filter(!is.na(region))
comb$point_fill <- ifelse(comb$sig_fdr,
                          PREDICTOR_COLOURS[as.character(comb$predictor)], "white")

# Delta-AIC per response x gene x region. Anchored to each panel's right edge
# because facet_grid uses free x-scales here, so a single computed x would land
# off-panel in some facets.
comb_lab <- comb |>
  distinct(resp_lab, subtype, region, AIC, AIC_null) |>
  mutate(txt = sprintf("ΔAIC = %.1f", AIC_null - AIC))

p_comb <- ggplot(comb, aes(x = estimate, y = region,
                           colour = predictor, group = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.6,
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = point_fill), shape = 21, size = 2.2, stroke = 0.7,
             position = position_dodge(width = 0.7)) +
  geom_text(data = comb_lab, aes(x = Inf, y = region, label = txt),
            hjust = 1.05, inherit.aes = FALSE, size = 2.1, colour = "grey40") +
  scale_colour_manual(values = PREDICTOR_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.34))) +
  facet_grid(resp_lab ~ subtype, scales = "free_x") +
  labs(title = "SLAC dN, dS and omega vs climate",
       subtitle = "Coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05",
       x = "Regression coefficient", y = NULL) +
  theme_bw(base_size = 11) +
  theme(strip.text = element_text(face = "bold", size = 9),
        strip.background = element_blank(),
        legend.position = "bottom",
        panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank()) +
  guides(colour = guide_legend(nrow = 2))

ggsave(file.path(fig_dir, "phylolm_slac_forest_combined.png"), p_comb,
       width = 9, height = 9, dpi = 200)
message("  Saved: phylolm_slac_forest_combined.png")

message("\nDone.")
