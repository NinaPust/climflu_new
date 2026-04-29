# scripts/03_rates_analysis/08_inla_country.R
#
# PURPOSE: INLA sensitivity analysis — restrict to countries with dense
#          year-round surveillance and add country + year fixed effects.
#          Mirrors 06_phylolm_country.R but uses INLA for full posteriors.
#
# COUNTRIES: USA, United Kingdom, Australia, China, Japan.
#            Edit DENSE_COUNTRIES below to adjust.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_inla_country_results.csv
#   data/rates/all_inla_country_results.csv
#
# FLAGS:
#   --gene X   run a single gene segment
#   --test     use subsampled data from data/test/

suppressPackageStartupMessages({
  if (!requireNamespace("INLA", quietly = TRUE))
    stop(
      "INLA is not installed. Run this in R first:\n\n",
      "  install.packages(\"INLA\",\n",
      "    repos = c(getOption(\"repos\"),\n",
      "              INLA = \"https://inla.r-inla-download.org/R/stable\"),\n",
      "    dep = TRUE)\n"
    )
  library(INLA)
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
  library(Matrix)
})

source(here("config.R"))
source(here("R", "climate_utils.R"))

# ── Settings ───────────────────────────────────────────────────────────────────
N_SUB <- 2000

DENSE_COUNTRIES <- c("USA", "United Kingdom", "Australia", "China", "Japan")

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args      <- commandArgs(trailingOnly = TRUE)
genes     <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}
use_test <- "--test" %in% args

if (use_test) message("*** Running in TEST mode ***")

paths_timed <- if (use_test) PATHS$test_trees_timed else PATHS$trees_timed
paths_rates <- if (use_test) PATHS$test_rates       else PATHS$rates

PREDICTOR_TERMS <- c("scaled_temp", "scaled_humidity",
                     "scaled_delta_temp", "scaled_delta_humidity",
                     "scaled_sampling_density", "scaled_pop_density")
                     # "scaled_season"

# ── Helpers (copied from 04_inla_phylo.R) ─────────────────────────────────────
posterior_summary <- function(fit, term) {
  marg <- fit$marginals.fixed[[term]]
  if (is.null(marg)) return(NULL)

  q <- inla.qmarginal(c(0.025, 0.975), marg)
  tibble(
    mean   = inla.emarginal(function(x) x, marg),
    sd     = sqrt(inla.emarginal(function(x) x^2, marg) -
                    inla.emarginal(function(x) x, marg)^2),
    q0.025 = q[1],
    q0.975 = q[2],
    p_pos  = 1 - inla.pmarginal(0, marg),
    sig_95 = as.integer(q[1] > 0 | q[2] < 0)
  )
}

stratified_subsample <- function(df, n, seed = 42) {
  set.seed(seed)
  if (nrow(df) <= n) return(df)
  total_n    <- nrow(df)
  year_counts <- df |>
    count(year, name = "grp_n") |>
    mutate(.k = pmax(1L, as.integer(round(grp_n / total_n * n))))
  sampled <- lapply(split(df, df$year), function(sub) {
    k <- year_counts$.k[year_counts$year == sub$year[1]]
    sub[sample(nrow(sub), min(k, nrow(sub))), ]
  })
  result <- do.call(rbind, sampled)
  result[sample(nrow(result), min(n, nrow(result))), ]
}

# ── Main loop ─────────────────────────────────────────────────────────────────
all_results <- list()

for (gene in genes) {
  message("\n[", gene, "] INLA (dense-country subset)...")

  rates_path <- file.path(paths_rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(paths_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(paths_rates,
                          paste0(gene, "_inla_country_results.csv"))

  if (!file.exists(rates_path)) { warning("Missing: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Missing: ", tree_path);  next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  # ── Restrict to dense-surveillance countries ────────────────────────────
  n_all <- nrow(df)
  df    <- df |> filter(country %in% DENSE_COUNTRIES)
  message("  Country filter: ", n_all, " → ", nrow(df), " rows (",
          paste(sort(unique(df$country)), collapse = ", "), ")")

  if (nrow(df) < 10) {
    warning("Too few rows after country filter for ", gene, ", skipping.")
    next
  }

  # ── Standard alignment ──────────────────────────────────────────────────
  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |>
    filter(accession %in% common, is.finite(log_rate)) |>
    arrange(match(accession, tree$tip.label))

  has_pop <- "scaled_pop_density" %in% names(df)
  if (has_pop) df <- df |> filter(!is.na(scaled_pop_density))

  # Drop sparse country/year levels
  df <- df |>
    group_by(country) |> filter(n() >= 2) |> ungroup() |>
    group_by(year)    |> filter(n() >= 2) |> ungroup()

  tree <- keep.tip(tree, df$accession)

  predictors <- c("scaled_temp", "scaled_humidity",
                  if ("scaled_delta_temp"       %in% names(df)) "scaled_delta_temp",
                  if ("scaled_delta_humidity"   %in% names(df)) "scaled_delta_humidity",
                  if ("scaled_sampling_density" %in% names(df)) "scaled_sampling_density",
                  if (has_pop) "scaled_pop_density")
                  # if ("scaled_season" %in% names(df)) "scaled_season"

  message("  ", nrow(df), " taxa available; subsampling to ",
          min(N_SUB, nrow(df)))

  # ── Subsample ───────────────────────────────────────────────────────────
  sub_df   <- stratified_subsample(df, N_SUB)
  sub_tree <- keep.tip(tree, sub_df$accession)
  sub_df   <- sub_df |>
    arrange(match(accession, sub_tree$tip.label))

  # ── Phylogenetic precision matrix ────────────────────────────────────────
  V <- vcv(sub_tree)
  V <- V[sub_df$accession, sub_df$accession]
  V <- V / max(diag(V))
  P <- Matrix(solve(V), sparse = FALSE)

  sub_df$tip_id_phylo <- seq_len(nrow(sub_df))
  sub_df$tip_id_iid   <- seq_len(nrow(sub_df))

  gene_results <- list()

  for (response in c("log_rate", "log_omega")) {
    if (!response %in% names(sub_df) ||
        all(is.na(sub_df[[response]]))) next

    year_term    <- if ("year"    %in% names(sub_df)) "as.factor(year)"    else NULL
    country_term <- if ("country" %in% names(sub_df)) "as.factor(country)" else NULL

    formula_inla <- as.formula(paste(
      response, "~",
      paste(c(predictors, year_term, country_term), collapse = " + "),
      "+ f(tip_id_phylo, model = 'generic0', Cmatrix = P,",
      "    hyper = list(prec = list(prior = 'pc.prec',",
      "                             param = c(1, 0.01))))",
      "+ f(tip_id_iid, model = 'iid',",
      "    hyper = list(prec = list(prior = 'pc.prec',",
      "                             param = c(1, 0.01))))"
    ))

    message("  Fitting INLA: ", response, " ~ ", gene)

    fit <- tryCatch(
      inla(
        formula_inla,
        data    = sub_df,
        family  = "gaussian",
        control.compute = list(dic = TRUE, waic = TRUE, cpo = FALSE),
        control.predictor = list(compute = FALSE),
        verbose = FALSE
      ),
      error = function(e) {
        warning("INLA failed for ", gene, " / ", response, ": ", e$message)
        NULL
      }
    )

    if (is.null(fit)) next

    # Extract only the focal predictors (year/country dummies are excluded)
    all_terms <- c("(Intercept)", predictors)
    rows <- lapply(all_terms, function(term) {
      s <- posterior_summary(fit, term)
      if (is.null(s)) return(NULL)
      bind_cols(
        tibble(gene = gene, response = response, term = term,
               n = nrow(sub_df)),
        s,
        tibble(dic = fit$dic$dic, waic = fit$waic$waic)
      )
    })

    gene_results[[response]] <- bind_rows(rows)
    message("  Done. DIC = ", round(fit$dic$dic, 1),
            "  WAIC = ", round(fit$waic$waic, 1))
  }

  result_df <- bind_rows(gene_results)
  write_csv(result_df, out_path)
  message("  Written: ", out_path)

  result_df |>
    filter(term %in% PREDICTOR_TERMS) |>
    mutate(across(c(mean, sd, q0.025, q0.975), \(x) round(x, 4)),
           p_pos = round(p_pos, 3)) |>
    print()

  all_results[[gene]] <- result_df
}

# ── Combine and summarise ─────────────────────────────────────────────────────
if (length(all_results) > 0) {
  combined      <- bind_rows(all_results)
  combined_path <- file.path(paths_rates, "all_inla_country_results.csv")
  write_csv(combined, combined_path)
  message("\nWritten combined: ", combined_path)

  message("\n── Predictor posteriors (log_rate, all genes) ──────────────────")
  combined |>
    filter(term %in% PREDICTOR_TERMS, response == "log_rate") |>
    mutate(
      ci    = paste0("[", round(q0.025, 3), ", ", round(q0.975, 3), "]"),
      p_pos = round(p_pos, 3),
      sig   = ifelse(sig_95 == 1, "*", "")
    ) |>
    select(gene, term, mean, ci, p_pos, sig) |>
    arrange(gene, term) |>
    print(n = 50)
}

message("\nDone.")
