# R/hurdle_utils.R
# Shared machinery for the two-part (hurdle) substitution analysis in
# scripts/04_epitope_analysis/12a, 12b, 12c.
#
# The hurdle splits one question into two that can be answered without dividing
# counts by branch time:
#
#   Part 1 (12a)  all tips            did ANY substitution occur?   -> rate
#   Part 2 (12b)  tips with 1 sub     was it protein-changing?      -> selection
#
# Part 1 needs an exposure covariate because a longer branch is mechanically
# more likely to carry a change. Part 2 does not: both outcomes required the
# same branch to exist, so branch time is common to numerator and denominator
# and cancels. That cancellation is the reason the hurdle avoids the artefact
# in which log(count / branch_time) is driven almost entirely by -log(time).

#' Cap BLAS/OpenMP threads.
#'
#' R here links against OpenBLAS (pthread build), which by default spawns one
#' thread per core. On a shared machine a single phyloglm fit was observed at
#' ~4100% CPU for no speed benefit: these fits are small and largely sequential,
#' so the extra threads are pure contention.
#'
#' Sys.setenv() alone is unreliable because OpenBLAS reads the variable when the
#' library loads, which happens before script code runs. RhpcBLASctl sets it
#' through the runtime API and works mid-session.
#'
#' @param n Maximum threads (default 4, or CLIMFLU_BLAS_THREADS).
#' @return The thread count applied, invisibly.
cap_blas_threads <- function(n = as.integer(Sys.getenv("CLIMFLU_BLAS_THREADS", "4"))) {
  if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
    RhpcBLASctl::blas_set_num_threads(n)
    RhpcBLASctl::omp_set_num_threads(n)
    message("BLAS/OMP threads capped at ", n, " (RhpcBLASctl)")
  } else if (!nzchar(Sys.getenv("OPENBLAS_NUM_THREADS"))) {
    warning("RhpcBLASctl not installed and OPENBLAS_NUM_THREADS unset — ",
            "OpenBLAS will use every core on this shared machine.\n",
            "  install.packages(\"RhpcBLASctl\")  or re-launch as:\n",
            "  OPENBLAS_NUM_THREADS=", n, " OMP_NUM_THREADS=", n, " Rscript <script>",
            immediate. = TRUE)
  } else {
    message("BLAS threads from environment: OPENBLAS_NUM_THREADS=",
            Sys.getenv("OPENBLAS_NUM_THREADS"))
  }
  invisible(n)
}


#' Terminal branch durations, measured on the UNPRUNED tree.
#'
#' Must be taken before any keep.tip()/drop.tip(), because pruning collapses
#' internal nodes and reassigns terminal branch lengths. Measuring after
#' pruning would make exposure incomparable between regions that happen to
#' retain different tip sets.
#'
#' @param tree phylo object (the full timetree).
#' @return Named numeric vector, branch duration per tip label.
terminal_branch_times <- function(tree) {
  idx <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
  stats::setNames(tree$edge.length[idx], tree$tip.label)
}


#' Align a data frame to a tree, dropping zero-length terminal branches.
#'
#' phyloglm requires one row per tip, in tip order, with row names matching
#' tip labels. Zero-length terminal branches are dropped because they make the
#' phylogenetic covariance singular.
#'
#' @param df    Data frame with an `accession` column.
#' @param tree  phylo object.
#' @param quiet Suppress the dropped-tip message.
#' @return list(data = <data.frame>, tree = <phylo>).
prep_tree_data <- function(df, tree, quiet = FALSE) {
  common <- intersect(tree$tip.label, df$accession)
  if (length(common) < 10)
    stop("Fewer than 10 tips shared between tree and data")
  tree <- ape::keep.tip(tree, common)

  idx       <- match(seq_len(ape::Ntip(tree)), tree$edge[, 2])
  zero_tips <- tree$tip.label[tree$edge.length[idx] == 0]
  if (length(zero_tips) > 0) {
    if (!quiet) message("    dropping ", length(zero_tips), " zero-length tips")
    tree <- ape::drop.tip(tree, zero_tips)
  }

  df <- df[df$accession %in% tree$tip.label, , drop = FALSE]
  df <- df[match(tree$tip.label, df$accession), , drop = FALSE]
  df <- as.data.frame(df)
  rownames(df) <- df$accession
  list(data = df, tree = tree)
}


#' Run a phyloglm Poisson-GEE count model and return a tidy coefficient table.
#'
#' Schema matches run_phyloglm() in phylolm_utils.R so results can be stacked
#' and plotted with shared code.
#'
#' GEE is a marginal (population-averaged) estimator: it uses the phylogeny to
#' build a working correlation structure and correct the standard errors,
#' rather than estimating a phylogenetic signal parameter. So there is no alpha
#' or lambda to report, and no likelihood — AIC and logLik come back NA. In
#' exchange it is robust to moderate misspecification of that structure, which
#' is what makes it a useful third leg rather than a competitor to the two
#' phyloglm fits.
#'
#' phyloglm has no offset() support (confirmed: no offset handling anywhere in
#' its source), so log-exposure enters as an ordinary covariate. Under a true
#' Poisson process its coefficient should land near 1; treat a value far from 1
#' as a diagnostic about the process, not a nuisance.
#'
#' @param formula  count_response ~ predictors (exposure included as a term).
#' @param data     Data frame, row names matching tree tip labels.
#' @param tree     phylo object.
#' @param subset   Label for this data subset.
#' @param response Label for the response variable.
#' @param gene     Gene segment label.
#' @return Tidy data frame, one row per coefficient.
run_phyloglm_gee <- function(formula, data, tree,
                             subset = "full", response = NA_character_,
                             gene = NA_character_) {
  tree <- ape::unroot(tree)
  fit  <- phylolm::phyloglm(formula, data = data, phy = tree,
                            method = "poisson_GEE")

  # Non-convergence must not pass silently. On sparse responses (e.g. the ridge,
  # ~2% of tips with a change) the GEE returns coefficients in the hundreds with
  # NaN standard errors, which would otherwise be written to the results CSV and
  # plotted as if they were estimates.
  if (!is.null(fit$convergence) && !identical(as.integer(fit$convergence), 0L))
    stop("poisson_GEE did not converge (convergence code ",
         fit$convergence, ") — response is probably too sparse for this region")

  coef_df <- as.data.frame(summary(fit)$coefficients)
  if (any(!is.finite(coef_df[[2]])))
    stop("poisson_GEE returned non-finite standard errors — ",
         "response is probably too sparse for this region")
  coef_df$term <- rownames(coef_df)
  rownames(coef_df) <- NULL

  # GEE summaries carry 4 columns (estimate, se, z, p); guard in case a
  # phylolm version returns fewer.
  if (ncol(coef_df) >= 5) {
    names(coef_df)[1:4] <- c("estimate", "std_error", "z_value", "p.value")
  } else {
    stop("Unexpected poisson_GEE summary layout: ", ncol(coef_df), " columns")
  }

  response_var <- all.vars(formula)[1]

  coef_df$gene      <- gene
  coef_df$subset    <- subset
  coef_df$response  <- response
  coef_df$alpha     <- NA_real_
  coef_df$logLik    <- NA_real_   # GEE is not likelihood-based
  coef_df$AIC       <- NA_real_
  coef_df$AIC_null  <- NA_real_
  coef_df$n         <- nrow(data)
  coef_df$n_events  <- sum(data[[response_var]] > 0, na.rm = TRUE)
  coef_df$r.squared <- NA_real_

  coef_df[, c("gene", "subset", "response", "term",
              "estimate", "std_error", "p.value",
              "alpha", "logLik", "AIC", "AIC_null", "n", "n_events", "r.squared")]
}


#' Apply BH-FDR correction across climate predictor rows only.
#'
#' Exposure, intercept and year dummies are excluded: they are nuisance terms,
#' not hypotheses, and including them would dilute the correction.
#'
#' @param res   Combined results data frame with `term` and `p.value`.
#' @param terms Character vector of predictor terms to correct over.
#' @return `res` with a `p.fdr` column (NA outside `terms`).
add_fdr <- function(res, terms) {
  res$p.fdr <- NA_real_
  rows <- res$term %in% terms
  if (any(rows)) res$p.fdr[rows] <- stats::p.adjust(res$p.value[rows], method = "BH")
  res
}
