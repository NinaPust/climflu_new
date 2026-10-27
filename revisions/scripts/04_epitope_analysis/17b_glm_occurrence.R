# revisions/scripts/04_epitope_analysis/17b_glm_occurrence.R
#
# PURPOSE: The same occurrence models as 17_phyloglm_slac_occurrence.R, fitted
#          as ordinary binomial GLMs instead of phyloglm.
#
# WHY: in 17, phyloglm's alpha settled on its UPPER bound in every non-ridge
#   model (alphaWarn = 2; e.g. H1 whole gene any_S reached 3.5341 against a
#   ceiling of 3.5349). alpha is the Ornstein-Uhlenbeck reversion rate, so the
#   upper bound means the fit wants alpha -> infinity, i.e. phylogenetic
#   correlation -> ZERO. phyloglm says as much in its own warning: "may simply
#   reflect a flat likelihood at large alpha values, meaning that the
#   phylogenetic correlation is estimated to be negligible."
#
#   On that flat surface the optimiser stops arbitrarily, and the full and null
#   fits settle in different places. The result was delta AIC values as low as
#   -519, which is impossible for nested models: with k parameters dropped,
#   delta AIC >= -2k, so the floor here is -12. Eleven of twenty models in 17
#   fell below it. Raising log.alpha.bound does not rescue them — at 6 the
#   likelihood stops being finite and the fit fails outright.
#
#   So the boundary behaviour is not a nuisance to tune away; it is the model
#   reporting that there is no phylogenetic signal left to correct for. Fitting
#   these as plain GLMs estimates the same quantity without the pathology, and
#   gives AICs that behave.
#
#   The ridge models are the exception and should stay with phyloglm: there
#   alpha was genuinely estimated (0.22-0.44, alphaWarn = 0) and the delta AICs
#   sat sensibly at -12.9, -12.7 and -4.6.
#
# THE HONEST LIMITATION: a GLM treats tips as independent, which they are not.
#   Point estimates are unaffected but standard errors will be too narrow. Two
#   things are therefore reported alongside:
#     - cluster-robust standard errors, clustering tips on their parent node
#       (CR1 sandwich, implemented here to avoid a dependency on `sandwich`).
#       This absorbs the correlation among immediate relatives, which is where
#       most of it sits once the phylogenetic term has been judged negligible.
#     - the dispersion statistic, as a check on the binomial assumption.
#
# MODEL (identical to 17):
#   any_N ~ 6 climate + log(EN * t) + as.factor(year)
#   any_S ~ 6 climate + log(ES * t) + as.factor(year)
#   null:  ~ log(exposure) + as.factor(year)     -> delta AIC isolates climate
#
# INPUT:
#   revisions/data/rates/{gene}_rates_slac.csv    from 15_slac_dnds.R
#   data/trees/timed/{gene}_timetree.nwk          durations and parent lookup
#
# OUTPUT:
#   revisions/data/rates/phylolm/{gene}_glm_occurrence_results.csv
#   revisions/data/rates/phylolm/all_glm_occurrence_results.csv
#   revisions/figures/epitope/glm_occurrence_forest{,_fdr}.png
#
# FLAGS:
#   --gene X        default: H1 H3
#   --region X      default: full HA1 HA2 epitope ridge
#   --response X    any_N | any_S (default: both)
#   --no-year       drop as.factor(year)
#   --no-robust     skip the cluster-robust standard errors
#   --plot-only     rebuild figures from existing per-gene CSVs
#
# RUNTIME: seconds. Each fit is a GLM on at most ~15k rows.

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
  library(tibble)
  library(ggplot2)
  library(viridisLite)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "hurdle_utils.R"))
source(here("revisions", "R", "forest_utils.R"))

args      <- commandArgs(trailingOnly = TRUE)
genes     <- if ("--gene"     %in% args) args[which(args == "--gene")     + 1] else c("H1", "H3")
regions   <- if ("--region"   %in% args) args[which(args == "--region")   + 1] else
             c("full", "HA1", "HA2", "epitope", "ridge")
responses <- if ("--response" %in% args) args[which(args == "--response") + 1] else
             c("any_N", "any_S")
use_year  <- !("--no-year"   %in% args)
do_robust <- !("--no-robust" %in% args)
plot_only <- "--plot-only"   %in% args

PRED <- c("scaled_temp", "scaled_humidity",
          "scaled_delta_temp", "scaled_delta_humidity",
          "scaled_sampling_density", "scaled_pop_density")

RESPONSE_LABELS <- c(any_N = "Occurrence of a non-synonymous change",
                     any_S = "Occurrence of a synonymous change")
RESPONSE_SPEC <- list(any_N = list(count = "N", sites = "EN"),
                      any_S = list(count = "S", sites = "ES"))

out_dir <- PATHS$rates_phylolm
fig_dir <- here("revisions", "figures", "epitope")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

gene_file <- function(g)
  file.path(out_dir, paste0(g, "_glm_occurrence_results.csv"))


#' Cluster-robust (CR1) standard errors for a GLM.
#'
#' V = B M B, with bread B = vcov(fit) and meat M the sum over clusters of the
#' outer product of summed score contributions. For a canonical link the score
#' for observation i is (y_i - mu_i) * x_i. The CR1 small-sample correction
#' G/(G-1) * (n-1)/(n-p) matches what sandwich::vcovCL applies by default.
#'
#' Implemented directly because `sandwich` is not installed here.
#'
#' @param fit     A fitted glm.
#' @param cluster Cluster id per row, aligned with the rows used in the fit.
#' @return Named vector of robust standard errors.
cluster_robust_se <- function(fit, cluster) {
  X  <- stats::model.matrix(fit)
  u  <- stats::residuals(fit, type = "response")
  if (length(cluster) != nrow(X))
    stop("cluster vector (", length(cluster), ") does not align with the ",
         "model matrix (", nrow(X), " rows)")
  sc <- X * u
  br <- stats::vcov(fit)

  idx  <- split(seq_len(nrow(sc)), cluster)
  meat <- matrix(0, ncol(X), ncol(X))
  for (i in idx) meat <- meat + tcrossprod(colSums(sc[i, , drop = FALSE]))

  G <- length(idx); n <- nrow(X); p <- ncol(X)
  adj <- (G / (G - 1)) * ((n - 1) / (n - p))
  sqrt(diag(br %*% (adj * meat) %*% br))
}

# ── Fit ───────────────────────────────────────────────────────────────────────

if (!plot_only) for (gene in genes) {
  message("\n========== ", gene, " ==========")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_slac.csv"))
  tree_path  <- here("data", "trees", "timed", paste0(gene, "_timetree.nwk"))
  if (!file.exists(rates_path)) {
    warning(gene, ": ", basename(rates_path), " not found — run 15_slac_dnds.R")
    next
  }
  if (!file.exists(tree_path)) { warning(gene, ": timetree not found"); next }

  df <- suppressWarnings(read_csv(rates_path, show_col_types = FALSE,
                                  guess_max = 50000))
  tr <- read.tree(tree_path)

  # Durations from the unpruned tree, and each tip's parent node as the
  # clustering unit for the robust standard errors.
  tip_edge <- match(seq_len(Ntip(tr)), tr$edge[, 2])
  df$t      <- stats::setNames(tr$edge.length[tip_edge], tr$tip.label)[df$accession]
  df$parent <- stats::setNames(tr$edge[tip_edge, 1], tr$tip.label)[df$accession]

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

      keep <- intersect(c("response", "log_exposure", "parent", PRED,
                          if (use_year) "year"), names(sub))
      sub  <- sub[, keep, drop = FALSE]
      sub  <- sub[stats::complete.cases(sub) & is.finite(sub$log_exposure), ,
                  drop = FALSE]

      n_ev <- sum(sub$response == 1)
      if (nrow(sub) < 100 || n_ev < 20 || n_ev > nrow(sub) - 20) {
        message(sprintf("  [%-7s %-5s] n = %5d  events = %4d — too few, skipping",
                        r, resp, nrow(sub), n_ev))
        next
      }

      preds     <- intersect(PRED, names(sub))
      year_term <- if (use_year && "year" %in% names(sub)) "as.factor(year)" else NULL
      rhs_full  <- paste(c(preds, "log_exposure", year_term), collapse = " + ")
      rhs_null  <- paste(c("log_exposure", year_term), collapse = " + ")

      fit  <- glm(as.formula(paste("response ~", rhs_full)),
                  data = sub, family = binomial)
      null <- glm(as.formula(paste("response ~", rhs_null)),
                  data = sub, family = binomial)

      cm <- summary(fit)$coefficients
      keep_terms <- intersect(c(PRED, "log_exposure"), rownames(cm))

      rse <- rep(NA_real_, length(keep_terms))
      if (do_robust) {
        rse_all <- tryCatch(cluster_robust_se(fit, sub$parent),
                            error = function(e) {
                              message("      robust SE failed: ",
                                      conditionMessage(e)); NULL })
        if (!is.null(rse_all)) rse <- unname(rse_all[keep_terms])
      }

      disp <- sum(residuals(fit, "pearson")^2) / fit$df.residual
      dAIC <- AIC(null) - AIC(fit)

      gene_rows[[paste(r, resp)]] <- tibble(
        gene = gene, subset = r, response = resp, term = keep_terms,
        estimate   = cm[keep_terms, "Estimate"],
        std_error  = cm[keep_terms, "Std. Error"],
        robust_se  = rse,
        z_value    = cm[keep_terms, "z value"],
        p.value    = cm[keep_terms, "Pr(>|z|)"],
        p_robust   = 2 * stats::pnorm(-abs(cm[keep_terms, "Estimate"] / rse)),
        AIC = AIC(fit), AIC_null = AIC(null),
        n = nrow(sub), n_events = n_ev,
        n_clusters = dplyr::n_distinct(sub$parent),
        dispersion = disp
      )

      tmp <- gene_rows[[paste(r, resp)]] |> filter(term == "scaled_temp")
      message(sprintf("  [%-7s %-5s] n = %5d  events = %4d  temp = %+.4f (p = %-9.3g robust p = %-9.3g)  dAIC = %+7.1f  disp = %.2f",
                      r, resp, nrow(sub), n_ev, tmp$estimate, tmp$p.value,
                      tmp$p_robust, dAIC, disp))
    }
  }

  if (length(gene_rows) == 0) { message("  no successful fits"); next }
  write_csv(bind_rows(gene_rows), gene_file(gene))
  message("  Written: ", basename(gene_file(gene)))
}

# ── Combine and FDR-correct ───────────────────────────────────────────────────

all_res <- lapply(genes, function(g) {
  p <- gene_file(g)
  if (file.exists(p)) suppressWarnings(read_csv(p, show_col_types = FALSE)) else NULL
}) |> bind_rows()

if (nrow(all_res) == 0) stop("No results. Run without --plot-only first.")

all_res$p.fdr <- NA_real_
pr <- all_res$term %in% PRED
all_res$p.fdr[pr] <- p.adjust(all_res$p.value[pr], method = "BH")
all_res$p_robust_fdr <- NA_real_
all_res$p_robust_fdr[pr] <- p.adjust(all_res$p_robust[pr], method = "BH")

combined <- file.path(out_dir, "all_glm_occurrence_results.csv")
write_csv(all_res, combined)
message("\nWritten: ", basename(combined))

# ── Report ────────────────────────────────────────────────────────────────────

message("\nDelta AIC per model (floor for 6 dropped terms is -12):")
all_res |>
  distinct(gene, subset, response, AIC, AIC_null, n, n_events, dispersion) |>
  mutate(dAIC = round(AIC_null - AIC, 1), dispersion = round(dispersion, 3)) |>
  select(gene, subset, response, n, n_events, dAIC, dispersion) |>
  arrange(dAIC) |>
  as.data.frame() |>
  print(row.names = FALSE)

message("\nTemperature coefficient (OR per ~10 C), naive vs cluster-robust:")
all_res |>
  filter(term == "scaled_temp") |>
  mutate(OR = round(exp(estimate), 3), estimate = round(estimate, 4),
         se = round(std_error, 4), rse = round(robust_se, 4),
         infl = round(robust_se / std_error, 2),
         p = signif(p.value, 3), p_rob = signif(p_robust, 3)) |>
  select(gene, subset, response, estimate, OR, se, rse, infl, p, p_rob) |>
  arrange(response, gene, subset) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ── Forest plots ──────────────────────────────────────────────────────────────

prep_plot <- function(res, se_col, p_col) {
  res |>
    filter(term %in% PRED) |>
    mutate(
      se_use    = .data[[se_col]],
      ci_lo     = estimate - 1.96 * se_use,
      ci_hi     = estimate + 1.96 * se_use,
      sig       = !is.na(.data[[p_col]]) & .data[[p_col]] < 0.05,
      predictor = factor(PREDICTOR_LABELS[term], levels = PREDICTOR_LABELS),
      region    = factor(REGION_LABELS[subset], levels = REGION_LEVELS),
      subtype   = ifelse(gene %in% c("H1", "N1"), "H1N1", "H3N2"),
      resp_lab  = factor(RESPONSE_LABELS[response], levels = RESPONSE_LABELS)
    ) |>
    filter(!is.na(region), is.finite(ci_lo), is.finite(ci_hi))
}

make_forest <- function(df, subtitle, file) {
  df$point_fill <- ifelse(df$sig,
                          PREDICTOR_COLOURS[as.character(df$predictor)], "white")
  daic <- df |>
    distinct(resp_lab, subtype, region, AIC, AIC_null) |>
    mutate(txt = sprintf("ΔAIC = %.1f", AIC_null - AIC))

  p <- ggplot(df, aes(x = estimate, y = region,
                      colour = predictor, group = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
               linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.6,
                   position = position_dodge(width = 0.7)) +
    geom_point(aes(fill = point_fill), shape = 21, size = 2.2, stroke = 0.7,
               position = position_dodge(width = 0.7)) +
    geom_text(data = daic, aes(x = Inf, y = region, label = txt),
              hjust = 1.05, inherit.aes = FALSE, size = 2.1, colour = "grey40") +
    scale_colour_manual(values = PREDICTOR_COLOURS, name = NULL) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev, drop = TRUE) +
    scale_x_continuous(expand = expansion(mult = c(0.06, 0.34))) +
    facet_grid(resp_lab ~ subtype, scales = "free_x") +
    labs(title = "Occurrence of substitution vs climate (binomial GLM)",
         subtitle = subtitle, x = "Log-odds coefficient", y = NULL) +
    theme_bw(base_size = 11) +
    theme(strip.text = element_text(face = "bold", size = 9),
          strip.background = element_blank(),
          legend.position = "bottom",
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank()) +
    guides(colour = guide_legend(nrow = 2))

  ggsave(file.path(fig_dir, file), p, width = 9, height = 7, dpi = 200)
  message("  Saved: ", file)
}

make_forest(prep_plot(all_res, "std_error", "p.fdr"),
            "Model-based SE  |  Filled = BH-FDR p < 0.05",
            "glm_occurrence_forest_fdr.png")

if (do_robust && any(is.finite(all_res$robust_se)))
  make_forest(prep_plot(all_res, "robust_se", "p_robust_fdr"),
              paste0("Cluster-robust SE (tips clustered on parent node)  |  ",
                     "Filled = BH-FDR p < 0.05"),
              "glm_occurrence_forest_robust.png")

message("\nDone.")
