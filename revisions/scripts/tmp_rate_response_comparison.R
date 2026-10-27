# revisions/scripts/tmp_rate_response_comparison.R
#
# TEMPORARY DIAGNOSTIC — not part of the pipeline.
#
# PURPOSE: Refit the published whole-gene phylolm against three different
#          constructions of "substitution rate", to see how much of the
#          published temperature coefficient depends on which one is used, and
#          produce a fig-1C-style forest plot for the cleanest of them.
#
#   A  treetime_rate   TreeTime's own per-branch relaxed-clock rate, from
#                      substitution_rates.tsv. A regularised estimate: branches
#                      with no mutations are shrunk toward the mean rather than
#                      being zero, so NO tips are dropped (~2.8% sit at the
#                      0.1x shrinkage floor). Only ~22% of its variance is
#                      branch time, versus ~52-57% for B and C.
#
#   B  pruned_over_time  terminal branch length of the pruned ML tree divided by
#                      the timetree branch length. Structurally the same as the
#                      published response, from a different tree: the pruned
#                      tree's terminal lengths are integer substitution counts
#                      / alignment length, so ~32% are exactly zero and drop out
#                      at the log.
#
#   C  published      log_rate as computed by 03_rates_analysis/01_compute_rates.R
#                      (divergence_tree branch / timetree branch). The submitted
#                      result, refit here for a like-for-like comparison.
#
# For each response the script reports the temperature coefficient both as
# published and with log(branch time) added as a covariate, plus how much of the
# response's variance log(branch time) alone explains. If a response is largely
# -log(branch time) in disguise, controlling for branch time moves the
# temperature coefficient substantially (B and C) or barely at all (A).
#
# FOREST PLOT: response A, "as published" form, styled to match panel C of
#          scripts/05_final_figs/fig_main.R — same predictors, colours, protein
#          x subtype layout, dodge, fill-if-FDR-significant convention and
#          delta-AIC annotations. One deliberate deviation: the delta-AIC label
#          x-position is derived from the data range rather than hardcoded at
#          0.05, because response A's coefficients are roughly half the
#          magnitude of the published ones and a fixed 0.05 falls off-panel.
#
# Model matches 03_rates_analysis/02_phylolm.R: the six scaled climate
# predictors + as.factor(year), phylolm with Pagel's lambda. The year-only null
# model (for delta AIC) is fitted only for response A's plotted variant.
#
# FDR: BH across all gene x predictor rows for response A. Note the published
#      figure's FDR was computed across genes AND responses (log_rate and
#      log_omega together), so p.fdr values are not directly comparable.
#
# INPUT (all from the ORIGINAL submitted tree, not revisions/):
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#   data/trees/pruned/{gene}_pruned.nwk
#   data/trees/timed/{gene}_treetime_out/substitution_rates.tsv
#
# OUTPUT:
#   revisions/data/rates/phylolm/tmp_rate_response_comparison.csv
#   revisions/data/rates/phylolm/tmp_responseA_coefficients.csv
#   revisions/figures/fig1C_responseA_forest.png / .pdf
#
# USAGE:
#   Rscript revisions/scripts/tmp_rate_response_comparison.R
#   Rscript revisions/scripts/tmp_rate_response_comparison.R --gene H1
#
# FLAGS:
#   --gene X          restrict to one gene (default: all six)
#   --plot-with-time  plot the "+ log(branch time)" variant instead
#   --skip-bc         fit response A only (skips the B/C diagnostic, faster)
#   --max-n N         subsample per fit (calibration only)
#   --threads N       BLAS/OMP cap (default 4)

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(phylolm)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "hurdle_utils.R"))

ALL_GENES <- c("H1", "H3", "N1", "N2", "PA_H1N1", "PA_H3N2")

args      <- commandArgs(trailingOnly = TRUE)
genes     <- if ("--gene"    %in% args) args[which(args == "--gene")    + 1] else ALL_GENES
plot_time <- "--plot-with-time" %in% args
skip_bc   <- "--skip-bc"     %in% args
max_n     <- if ("--max-n"   %in% args) as.integer(args[which(args == "--max-n")   + 1]) else NA_integer_
n_thread  <- if ("--threads" %in% args) as.integer(args[which(args == "--threads") + 1]) else 4L

cap_blas_threads(n_thread)
set.seed(42)

PRED <- c("scaled_temp", "scaled_humidity",
          "scaled_delta_temp", "scaled_delta_humidity",
          "scaled_sampling_density", "scaled_pop_density")

root <- here()

RESPONSES <- list(
  A = "A treetime per-branch rate (no zeros)",
  B = "B pruned tree / timetree",
  C = "C published log_rate (divergence/timetree)"
)
if (skip_bc) RESPONSES <- RESPONSES["A"]


#' Build the per-tip table of all three responses for one gene.
build_gene_data <- function(gene) {
  rates_csv <- file.path(root, "data", "rates",  paste0(gene, "_rates.csv"))
  tt_path   <- file.path(root, "data", "trees", "timed",  paste0(gene, "_timetree.nwk"))
  pruned_p  <- file.path(root, "data", "trees", "pruned", paste0(gene, "_pruned.nwk"))
  subrates  <- file.path(root, "data", "trees", "timed",
                         paste0(gene, "_treetime_out"), "substitution_rates.tsv")
  miss <- c(rates_csv, tt_path, pruned_p, subrates)
  miss <- miss[!file.exists(miss)]
  if (length(miss)) { warning(gene, ": missing ", paste(basename(miss), collapse=", ")); return(NULL) }

  df <- read_csv(rates_csv, show_col_types = FALSE)
  tt <- read.tree(tt_path)
  pr <- read.tree(pruned_p)

  sr <- read.delim(subrates, comment.char = "")
  names(sr)[1] <- "node"

  df$time_bl <- terminal_branch_times(tt)[df$accession]
  df$pr_bl   <- terminal_branch_times(pr)[df$accession]
  df$tt_rate <- sr$rate[match(df$accession, sr$node)]

  df$resp_A <- ifelse(!is.na(df$tt_rate) & df$tt_rate > 0, log(df$tt_rate), NA_real_)
  df$resp_B <- ifelse(!is.na(df$pr_bl) & df$pr_bl > 0 & df$time_bl > 0,
                      log(df$pr_bl / df$time_bl), NA_real_)
  df$resp_C <- df$log_rate
  df$logt   <- ifelse(!is.na(df$time_bl) & df$time_bl > 0, log(df$time_bl), NA_real_)

  list(df = df, tree = tt)
}


#' Fit one response for one gene. Returns the full tidy coefficient table.
fit_one <- function(gd, col, add_time, want_null = FALSE) {
  df <- gd$df
  keep <- c("accession", col, "logt", "year", PRED)
  d <- df[, keep, drop = FALSE]
  d <- d[stats::complete.cases(d) & is.finite(d[[col]]) & is.finite(d$logt), , drop = FALSE]
  if (nrow(d) < 50) return(NULL)
  if (!is.na(max_n) && nrow(d) > max_n) d <- d[sample(nrow(d), max_n), ]

  p  <- prep_tree_data(d, gd$tree, quiet = TRUE)
  d  <- p$data; tr <- ape::unroot(p$tree)

  terms <- c(PRED, if (add_time) "logt", "as.factor(year)")
  fml   <- as.formula(paste(col, "~", paste(terms, collapse = " + ")))

  t0  <- proc.time()
  fit <- tryCatch(phylolm(fml, data = d, phy = tr, model = PARAMS$phylolm_method),
                  error = function(e) { message("      FAILED: ", conditionMessage(e)); NULL })
  if (is.null(fit)) return(NULL)
  el <- (proc.time() - t0)[["elapsed"]]

  cf <- as.data.frame(summary(fit)$coefficients)
  names(cf)[1:4] <- c("estimate", "std_error", "t_value", "p.value")
  cf$term <- rownames(cf); rownames(cf) <- NULL

  aic_null <- NA_real_
  if (want_null) {
    nf <- as.formula(paste(col, "~ as.factor(year)"))
    aic_null <- tryCatch(AIC(phylolm(nf, data = d, phy = tr,
                                     model = PARAMS$phylolm_method)),
                         error = function(e) NA_real_)
  }

  cf$AIC      <- AIC(fit)
  cf$AIC_null <- aic_null
  cf$lambda   <- tryCatch(fit$optpar, error = function(e) NA_real_)
  cf$n        <- nrow(d)
  cf$elapsed  <- el
  # Variance of the response explained by branch time alone (plain lm).
  cf$r2_time  <- summary(lm(as.formula(paste(col, "~ logt")), data = d))$r.squared
  cf
}

# ── Main loop ─────────────────────────────────────────────────────────────────

summary_rows <- list()
respA_coefs  <- list()

for (gene in genes) {
  message("\n========== ", gene, " ==========")
  gd <- build_gene_data(gene)
  if (is.null(gd)) next
  message("  tips in rates CSV: ", nrow(gd$df),
          "  matched to substitution_rates.tsv: ", sum(!is.na(gd$df$tt_rate)))

  for (k in names(RESPONSES)) {
    col     <- paste0("resp_", k)
    n_avail <- sum(is.finite(gd$df[[col]]))
    message("\n  [", k, "] ", RESPONSES[[k]],
            "\n      finite: ", n_avail, " of ", nrow(gd$df),
            "  (dropped ", nrow(gd$df) - n_avail, ")")

    for (add_time in c(FALSE, TRUE)) {
      lab  <- if (add_time) "+ log(branch time)" else "as published"
      # Null fit only where it is consumed: response A's plotted variant.
      want_null <- (k == "A") && (add_time == plot_time)
      cf <- fit_one(gd, col, add_time, want_null = want_null)
      if (is.null(cf)) { message("      ", lab, ": no fit"); next }

      tr <- cf[cf$term == "scaled_temp", ]
      message(sprintf("      %-20s n=%5d  temp=%+.4f (se %.4f)  p=%-10.3g lambda=%.3f  [%.0fs]",
                      lab, tr$n, tr$estimate, tr$std_error, tr$p.value, tr$lambda, tr$elapsed))
      if (add_time) {
        lt <- cf[cf$term == "logt", ]
        message(sprintf("      %-20s log(time) coefficient = %+.3f", "", lt$estimate))
      } else {
        message(sprintf("      %-20s response variance explained by log(time) alone: %.1f%%",
                        "", 100 * tr$r2_time))
      }

      summary_rows[[paste(gene, k, add_time)]] <- data.frame(
        gene = gene, response = RESPONSES[[k]], model = lab, n = tr$n,
        temp = tr$estimate, se = tr$std_error, p = tr$p.value,
        logt_coef = if (add_time) cf$estimate[cf$term == "logt"] else NA_real_,
        lambda = tr$lambda, r2_time = tr$r2_time
      )

      if (k == "B" && add_time == plot_time) {
        cf$gene <- gene
        respA_coefs[[gene]] <- cf
      }
    }
  }
}

if (length(summary_rows) == 0) stop("No fits succeeded.")

res <- bind_rows(summary_rows)
out_dir <- PATHS$rates_phylolm
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
write_csv(res, file.path(out_dir, "tmp_rate_response_comparison.csv"))

message("\n================ SUMMARY (temperature coefficient) ================")
res |>
  mutate(temp = round(temp, 4), se = round(se, 4), p = signif(p, 3),
         logt_coef = round(logt_coef, 3), pct_var_time = round(100 * r2_time, 1)) |>
  select(gene, response, model, n, temp, se, p, logt_coef, pct_var_time) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ── Forest plot — response A, styled to match fig_main.R panel C ──────────────

if (length(respA_coefs) == 0) {
  message("\nNo response-A fits to plot.")
} else {

BS  <- 10   # base font size, as fig_main.R
BS2 <- 11

common_theme <- theme_bw(base_size = BS) +
  theme(
    strip.background   = element_blank(),
    strip.text         = element_text(face = "bold", size = BS + 1),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = BS2),
    axis.text          = element_text(size = BS - 1),
    legend.text        = element_text(size = BS - 1),
    legend.key.size    = unit(0.6, "lines"),
    plot.tag           = element_text(face = "bold.italic", size = BS + 2),
    plot.tag.position  = "topleft"
  )

KEEP_TERMS <- c("scaled_temp", "scaled_humidity",
                "scaled_sampling_density", "scaled_delta_temp",
                "scaled_delta_humidity", "scaled_pop_density")
TERM_LABELS <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative\nhumidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative\nhumidity",
  scaled_sampling_density = "Sampling\ndensity",
  scaled_pop_density      = "Population\ndensity"
)
PRED_COLOURS <- setNames(viridis(length(TERM_LABELS), end = 0.92), TERM_LABELS)

allA <- bind_rows(respA_coefs)
# BH-FDR across all gene x predictor rows (see header note on comparability).
allA$p.fdr <- NA_real_
pr_rows <- allA$term %in% KEEP_TERMS
allA$p.fdr[pr_rows] <- p.adjust(allA$p.value[pr_rows], method = "BH")
write_csv(allA, file.path(out_dir, "tmp_responseA_coefficients.csv"))

res_A <- allA |>
  filter(term %in% KEEP_TERMS) |>
  mutate(
    ci_lo   = estimate - 1.96 * std_error,
    ci_hi   = estimate + 1.96 * std_error,
    sig_fdr = !is.na(p.fdr) & p.fdr < 0.05,
    protein = factor(
      case_when(gene %in% c("H1", "H3") ~ "HA",
                gene %in% c("N1", "N2") ~ "NA",
                TRUE                    ~ "PA"),
      levels = c("HA", "NA", "PA")
    ),
    subtype   = ifelse(gene %in% c("H1", "N1", "PA_H1N1"), "H1N1", "H3N2"),
    predictor = factor(TERM_LABELS[term], levels = TERM_LABELS),
    pt_fill   = ifelse(sig_fdr, PRED_COLOURS[as.character(predictor)], "white")
  )

r2_forest <- res_A |>
  distinct(gene, subtype, protein, AIC, AIC_null) |>
  filter(!is.na(AIC_null)) |>
  mutate(r2_label = sprintf("ΔAIC = %.1f", AIC_null - AIC))

# Deviation from fig_main.R: label x fixed at 0.05 there; here derived from the
# data range, because response A's coefficients are ~half the published size.
lab_x <- max(res_A$ci_hi, na.rm = TRUE) * 0.55
r2_forest$lab_x <- lab_x

dodge         <- position_dodge(width = 0.5)
forest_labels <- data.frame(subtype = unique(res_A$subtype))

p_forest <- ggplot(res_A,
                   aes(x = estimate, y = protein,
                       colour = predictor, group = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 linewidth = 0.6, position = dodge) +
  geom_point(aes(fill = pt_fill),
             shape = 21, size = 2.2, stroke = 0.7, position = dodge) +
  geom_text(data = forest_labels,
            aes(x = -Inf, y = Inf, label = subtype),
            hjust = -0.1, vjust = 1.4, inherit.aes = FALSE,
            size = (BS - 1) / .pt, colour = "grey20", fontface = "bold") +
  geom_text(data = r2_forest,
            aes(x = lab_x, y = protein, label = r2_label),
            hjust = 0, inherit.aes = FALSE,
            size = (BS - 3) / .pt, colour = "grey40") +
  scale_colour_manual(values = PRED_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev) +
  facet_wrap(~ subtype, ncol = 1) +
  labs(x = "Regression coefficient", y = NULL, tag = "C") +
  common_theme +
  theme(legend.position      = c(-0.01, 0.01),
        legend.justification = c(0, 0),
        legend.direction     = "vertical",
        legend.background    = element_rect(fill = "transparent", colour = NA),
        legend.text          = element_text(size = BS - 3),
        legend.key.spacing.y = unit(0.1, "lines"),
        panel.grid.major.x   = element_blank(),
        aspect.ratio         = 0.5,
        strip.text           = element_blank(),
        strip.background     = element_blank()) +
  guides(colour = guide_legend(nrow = 6))

fig_dir <- here("revisions", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
stem <- if (plot_time) "fig1C_responseA_forest_withtime" else "fig1C_responseA_forest"
for (ext in c("png", "pdf"))
  ggsave(file.path(fig_dir, paste0(stem, ".", ext)), p_forest,
         width = 4.2, height = 4.6, dpi = 300)
message("\nSaved: ", stem, ".png / .pdf  in revisions/figures/")

message("\nResponse A, significant after BH-FDR:")
sigA <- res_A |>
  filter(sig_fdr) |>
  mutate(estimate = round(estimate, 4), p.fdr = signif(p.fdr, 3)) |>
  select(gene, protein, subtype, term, estimate, p.fdr, n) |>
  arrange(p.fdr)
if (nrow(sigA) == 0) message("  (none)") else print(as.data.frame(sigA), row.names = FALSE)

}

message("\nDone.")
