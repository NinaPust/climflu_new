# scripts/04_epitope_analysis/13_jackknife_ns_rates.R
#
# PURPOSE: Block jackknife subsampling for phylolm temperature analysis of
#          NS substitution rates across four HA regions (ridge, epitope,
#          full gene, HA2). Mirrors scripts/03_rates_analysis/07_phylolm_jackknife.R
#          but uses data/rates/{gene}_rates_epitope.csv and loops over regions.
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/jackknife/epitope/{gene}_{region}_jackknife_results.csv
#   data/rates/jackknife/epitope/{gene}_{region}_jackknife_membership.csv
#
# FLAGS:
#   --gene X     run a single gene (default: H1 H3)
#   --region X   run a single region (default: full HA2 epitope ridge)

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
  library(ape)
  library(dplyr)
  library(readr)
})

source(here("config.R"))

N_SAMPLES   <- 1000L
SAMPLE_SIZE <- 1000L

# Per-region seeds — offset from 42 so each region gives different samples
REGION_SEEDS <- c(full = 42L, HA2 = 43L, epitope = 44L, ridge = 45L)

args    <- commandArgs(trailingOnly = TRUE)
genes   <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
regions <- if ("--region" %in% args) args[which(args == "--region") + 1] else c("full", "HA2", "epitope", "ridge")

out_dir <- file.path(PATHS$rates, "jackknife", "epitope")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

for (gene in genes) {
  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))

  if (!file.exists(rates_path)) { warning("Rates file not found: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Tree not found: ",        tree_path);  next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  # Align tips and data; drop zero-length terminal branches
  common       <- intersect(tree$tip.label, df$accession)
  tree         <- keep.tip(tree, common)
  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    message("[", gene, "] Dropping ", length(zero_tips), " zero-length tips")
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }
  df <- df |> filter(accession %in% common)

  for (region in regions) {
    response <- paste0("log_", region, "_ns_rate")
    tag      <- paste0("[", gene, " ", region, "]")

    if (!response %in% names(df)) {
      warning(tag, " column '", response, "' not found in rates — skipping")
      next
    }

    message("\n", tag, " Jackknife (", N_SAMPLES, " × ", SAMPLE_SIZE, " sequences)...")

    sub_df <- df |>
      filter(!is.na(scaled_temp), !is.na(.data[[response]])) |>
      arrange(match(accession, tree$tip.label))
    sub_tree <- keep.tip(tree, sub_df$accession)
    rownames(sub_df) <- sub_df$accession

    n_avail     <- nrow(sub_df)
    actual_size <- min(SAMPLE_SIZE, n_avail)
    message("  ", n_avail, " taxa available; sampling ", actual_size)

    if (n_avail < 10) {
      warning(tag, " Too few taxa for jackknife, skipping")
      next
    }
    if (actual_size < SAMPLE_SIZE) {
      warning(tag, " Only ", n_avail, " sequences available; sampling all")
    }

    set.seed(REGION_SEEDS[[region]])

    results_list    <- vector("list", N_SAMPLES)
    membership_list <- vector("list", N_SAMPLES)

    for (i in seq_len(N_SAMPLES)) {
      sampled_acc      <- sample(sub_df$accession, actual_size, replace = FALSE)
      samp_df          <- sub_df |> filter(accession %in% sampled_acc)
      rownames(samp_df) <- samp_df$accession
      samp_tree        <- unroot(keep.tip(sub_tree, sampled_acc))

      formula_str  <- paste(response, "~ scaled_temp + as.factor(year)")

      fit <- tryCatch(
        phylolm(
          as.formula(formula_str),
          data  = samp_df,
          phy   = samp_tree,
          model = PARAMS$phylolm_method
        ),
        error = function(e) {
          message("  Sample ", i, " failed: ", conditionMessage(e))
          NULL
        }
      )

      if (is.null(fit)) next

      coef_mat <- summary(fit)$coefficients
      if (!"scaled_temp" %in% rownames(coef_mat)) next

      se_col <- intersect(c("StdErr", "Std.Error", "std_error"),
                          colnames(coef_mat))[1]

      results_list[[i]] <- tibble(
        gene      = gene,
        region    = region,
        sample_id = i,
        n         = actual_size,
        estimate  = coef_mat["scaled_temp", "Estimate"],
        std_error = coef_mat["scaled_temp", se_col],
        t_value   = coef_mat["scaled_temp", "t.value"],
        p.value   = coef_mat["scaled_temp", "p.value"],
        lambda    = fit$optpar
      )

      membership_list[[i]] <- tibble(
        sample_id = i,
        accession = sampled_acc
      )

      if (i %% 100 == 0) message("  Sample ", i, "/", N_SAMPLES, " done")
    }

    results_df    <- bind_rows(results_list)
    membership_df <- bind_rows(membership_list)

    res_path <- file.path(out_dir, paste0(gene, "_", region,
                                          "_jackknife_results.csv"))
    mem_path <- file.path(out_dir, paste0(gene, "_", region,
                                          "_jackknife_membership.csv"))
    write_csv(results_df,    res_path)
    write_csv(membership_df, mem_path)

    message("  Written: ", basename(res_path))
    message("  Estimate range: [",
            round(min(results_df$estimate, na.rm = TRUE), 4), ", ",
            round(max(results_df$estimate, na.rm = TRUE), 4), "]")
  }
}

message("\nDone.")
