# R/phylolm_utils.R
# Utilities for running phylolm models and formatting results.

#' Run a phylolm model and return a tidy coefficient table.
#'
#' @param formula  Model formula (response ~ predictors).
#' @param data     Data frame with row names matching tree tip labels.
#' @param tree     phylo object.
#' @param method   phylolm method (default "lambda").
#' @param subset   Label for this data subset (e.g. "full", "urban", "rural").
#' @param response Label for response variable (for output table).
#' @param gene     Gene segment label (for output table).
#' @return Tidy data frame with one row per coefficient.
run_phylolm <- function(formula, data, tree, method = "lambda",
                        subset = "full", response = NA_character_,
                        gene = NA_character_) {
  tree <- ape::unroot(tree)
  fit <- phylolm::phylolm(formula, data = data, phy = tree, model = method)

  coef_df <- as.data.frame(summary(fit)$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL

  names(coef_df) <- c("estimate", "std_error", "t_value", "p.value", "term")

  coef_df$gene      <- gene
  coef_df$subset    <- subset
  coef_df$response  <- response
  coef_df$lambda    <- tryCatch(fit$optpar, error = function(e) NA_real_)
  coef_df$logLik    <- as.numeric(logLik(fit))
  coef_df$AIC       <- AIC(fit)
  coef_df$n         <- nrow(data)
  coef_df$r.squared <- tryCatch(summary(fit)$r.squared, error = function(e) NA_real_)

  # Null model: year factor only (no climate predictors) — for delta AIC.
  response_var <- all.vars(formula)[1]
  has_year     <- "year" %in% names(data)
  null_formula <- if (has_year) {
    as.formula(paste(response_var, "~ as.factor(year)"))
  } else {
    as.formula(paste(response_var, "~ 1"))
  }
  coef_df$AIC_null <- tryCatch(
    AIC(phylolm::phylolm(null_formula, data = data, phy = tree, model = method)),
    error = function(e) NA_real_
  )

  coef_df[, c("gene", "subset", "response", "term",
              "estimate", "std_error", "t_value", "p.value",
              "lambda", "logLik", "AIC", "AIC_null", "n", "r.squared")]
}


#' Format a list of phylolm result tables into a publication-ready data frame.
#'
#' @param models Named list of data frames (output of run_phylolm).
#' @return Combined data frame, sorted by gene and subset.
format_phylolm_table <- function(models) {
  df <- dplyr::bind_rows(models)
  df <- dplyr::arrange(df, .data$gene, .data$subset, .data$response, .data$term)
  df
}


#' Run a phyloglm (binary presence/absence) model and return a tidy coefficient
#' table with the same column layout as run_phylolm(), so results from the two
#' methods can be combined/plotted with shared code.
#'
#' Response should be a 0/1 indicator (e.g. "any NS change in this region"),
#' not a rate. Uses the Ives & Garland (2010) penalized-likelihood method
#' ("logistic_MPLE"), which is more robust than plain logistic regression to
#' the rare-event / quasi-separation issues expected for sparse regions like
#' the antigenic ridge.
#'
#' @param formula  Model formula (binary_response ~ predictors).
#' @param data     Data frame with row names matching tree tip labels.
#' @param tree     phylo object.
#' @param btol     Bound on the boundary parameter passed to phyloglm (higher
#'                 values help convergence for rare/quasi-separated outcomes).
#' @param subset   Label for this data subset (e.g. "full").
#' @param response Label for response variable (for output table).
#' @param gene     Gene segment label (for output table).
#' @return Tidy data frame with one row per coefficient.
#' @param fit_null Fit a second, reduced model to obtain AIC_null for a
#'   delta-AIC comparison. This DOUBLES runtime, and phyloglm fits are the
#'   expensive step, so set FALSE unless the value is actually consumed.
#' @param null_terms Extra terms to KEEP in the null model alongside the year
#'   factor, as a character vector. Delta AIC is only interpretable as "what do
#'   the climate predictors add" if the null contains everything except them —
#'   so any nuisance covariate in the full model, in particular an exposure or
#'   offset-substitute term such as log(sites * branch time), must be named here.
#'   Omitting it makes the null a worse fit for reasons that have nothing to do
#'   with climate and inflates delta AIC accordingly. Default NULL reproduces the
#'   original year-only null.
#' @param log_alpha_bound Passed to phyloglm as `log.alpha.bound`. phyloglm
#'   bounds the phylogenetic correlation parameter at roughly
#'   exp(log.alpha.bound) / tree_height; its default of 4 puts the ceiling near
#'   3.5 for H1 and 5.5 for H3, which every non-ridge model in script 17 hit
#'   exactly. At the bound the likelihood surface is flat, so the optimiser
#'   stops arbitrarily and the full and null fits can settle in different
#'   places — which is how delta AIC reached -519, far below the -2k floor that
#'   nested maximum-likelihood models obey. Raise this when alpha pegs.
#'   The SAME bound is applied to the null fit; otherwise the two AICs are not
#'   comparable in the first place.
run_phyloglm <- function(formula, data, tree, btol = 30,
                         subset = "full", response = NA_character_,
                         gene = NA_character_, fit_null = TRUE,
                         null_terms = NULL, log_alpha_bound = 4) {
  tree <- ape::unroot(tree)
  fit <- phylolm::phyloglm(formula, data = data, phy = tree,
                           method = "logistic_MPLE", btol = btol,
                           log.alpha.bound = log_alpha_bound)

  s <- summary(fit)
  coef_df <- as.data.frame(s$coefficients)
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL

  names(coef_df) <- c("estimate", "std_error", "z_value", "p.value", "term")

  response_var <- all.vars(formula)[1]
  n_events     <- sum(data[[response_var]] == 1, na.rm = TRUE)

  coef_df$gene      <- gene
  coef_df$subset    <- subset
  coef_df$response  <- response
  coef_df$alpha     <- tryCatch(unname(fit$alpha), error = function(e) NA_real_)
  # phyloglm's own flag for whether alpha settled on a bound:
  #   0 = estimated freely, 1 = at the lower bound, 2 = at the upper bound.
  # A 1 or 2 means the likelihood was flat there and the fit is unreliable —
  # this is what produced delta AIC values below the nested-model floor.
  coef_df$alpha_warn <- tryCatch(unname(fit$alphaWarn), error = function(e) NA_integer_)
  # The bounds themselves are exp(+/- log.alpha.bound) / mean tip height.
  mth <- tryCatch(unname(fit$mean.tip.height), error = function(e) NA_real_)
  coef_df$alpha_upper <- exp(log_alpha_bound)  / mth
  coef_df$alpha_lower <- exp(-log_alpha_bound) / mth
  coef_df$logLik    <- as.numeric(logLik(fit))
  coef_df$AIC       <- AIC(fit)
  coef_df$n         <- nrow(data)
  coef_df$n_events  <- n_events
  coef_df$r.squared <- NA_real_  # not defined for phyloglm; kept for schema parity

  coef_df$AIC_null <- NA_real_
  if (isTRUE(fit_null)) {
    rhs <- c(null_terms,
             if ("year" %in% names(data)) "as.factor(year)")
    null_formula <- as.formula(
      paste(response_var, "~",
            if (length(rhs)) paste(rhs, collapse = " + ") else "1")
    )
    coef_df$AIC_null <- tryCatch(
      AIC(phylolm::phyloglm(null_formula, data = data, phy = tree,
                            method = "logistic_MPLE", btol = btol,
                            log.alpha.bound = log_alpha_bound)),
      error = function(e) NA_real_
    )
  }

  coef_df[, c("gene", "subset", "response", "term",
              "estimate", "std_error", "p.value",
              "alpha", "alpha_warn", "alpha_lower", "alpha_upper",
              "logLik", "AIC", "AIC_null", "n", "n_events", "r.squared")]
}
