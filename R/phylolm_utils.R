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
