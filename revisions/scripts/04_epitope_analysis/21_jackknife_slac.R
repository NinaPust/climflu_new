# revisions/scripts/04_epitope_analysis/21_jackknife_slac.R
#
# PURPOSE: Block-jackknife robustness check on the phylolm models in
#          16_phylolm_slac.R. Repeatedly refits the model on random subsamples
#          of the available branches and records the coefficients, so the
#          question "does this result depend on a handful of sequences?" has an
#          answer. Revised counterpart to 13_jackknife_ns_rates.R +
#          14_plot_jackknife_ns.R, merged into one script.
#
#          Same data and same formula as script 16 — SLAC-derived responses,
#          all six climate predictors plus the year factor. The original
#          script 13 used a reduced formula (temperature + year only); this one
#          does not, so the jackknifed coefficient means the same thing as the
#          one reported in figure 1D.
#
# NO FDR, NO NULL MODEL. Neither is meaningful here: the 1000 replicates are
#   deliberately overlapping subsamples of one dataset, not independent tests,
#   so correcting across them would be wrong; and delta AIC per replicate would
#   answer a question nobody asked. Significance colouring in the plots uses
#   raw p-values, purely as a visual guide to replicate-level stability.
#
# ALL COEFFICIENTS ARE STORED, not just the one being plotted. Results are saved
#   in long format (one row per replicate per term), so --coef can be changed
#   and the figures rebuilt with --plot-only, without refitting anything.
#
# CUSTOMISABLE ON THREE AXES:
#   --region    which part of the gene     (full | HA1 | HA2 | epitope | ridge)
#   --response  which mutation type        (log_dN_rate | log_dS_rate | log_omega)
#   --coef      which coefficient to plot  (any scaled_* predictor)
#
# RESPONSES (identical to 16_phylolm_slac.R):
#   log_dN_rate = log( N / (EN * t) )        requires N > 0
#   log_dS_rate = log( S / (ES * t) )        requires S > 0
#   log_omega   = log( (N/EN) / (S/ES) )     requires N > 0 AND S > 0; no t,
#                                            which cancels between the two
#
# A NOTE ON SAMPLE SIZE: the full formula costs about 17 parameters. With the
#   default 60% subsample capped at 1000 branches this is comfortable for the
#   whole gene and HA1/HA2, but the ridge has only ~125-187 branches with a
#   defined dN, so replicates there carry ~75-110 observations and will be
#   visibly unstable. That instability is the finding, not a failure — but it
#   is worth saying so explicitly rather than letting a wide spread look like
#   a result.
#
# INPUT:
#   revisions/data/rates/{gene}_rates_slac.csv    from 15_slac_dnds.R
#   data/trees/timed/{gene}_timetree.nwk          branch durations (years)
#
# OUTPUT:
#   revisions/data/rates/jackknife/slac/{draw}/{gene}_{region}_{response}_results.csv
#   revisions/data/rates/jackknife/slac/{draw}/..._membership.csv
#     where {draw} is n1000 / pct15 / pct30 / pct60. Unlike 03/07, which keeps
#     n1000 at its top level for backward compatibility, every draw size here
#     gets its own directory.
#   revisions/figures/jackknife/slac/{gene}_{region}_{response}_{coef}_samples.png
#   revisions/figures/jackknife/slac/..._samples_ordered.png
#   revisions/figures/jackknife/slac/..._sequences.png
#   revisions/figures/jackknife/slac/jackknife_summary_{response}_{coef}.png
#
# FLAGS:
#   --gene X        default: H1 H3
#   --region X      default: epitope
#   --response X    default: log_dN_rate
#   --coef X        coefficient to plot, default: scaled_temp
#   --samples N     jackknife replicates per draw size, default 1000
#   --draw X        which draw sizes to run: n1000 | pct15 | pct30 | pct60
#                   (default: all four, matching 03/07_phylolm_jackknife.R)
#   --plot-draw X   which draw size to plot, default pct60
#   --force         refit even if a results file already exists
#   --min-n N       skip a combination with fewer than N branches, default 40
#   --seed N        base seed, default 42
#   --threads N     BLAS/OMP cap, default 4
#   --plot-only     rebuild figures from existing results, no refitting
#
# RUNTIME: one phylolm fit per replicate per draw size, so roughly 4000 fits
#   per gene x region x response with the default four sizes. Results are
#   written per draw size and skipped if present, so an interrupted run resumes
#   where it stopped. Run it in tmux.

suppressPackageStartupMessages({
  library(here)
  library(phylolm)
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
regions   <- if ("--region"   %in% args) args[which(args == "--region")   + 1] else "epitope"
responses <- if ("--response" %in% args) args[which(args == "--response") + 1] else "log_dN_rate"
coef_plot <- if ("--coef"     %in% args) args[which(args == "--coef")     + 1] else "scaled_temp"
n_samples <- if ("--samples"  %in% args) as.integer(args[which(args == "--samples")  + 1]) else 1000L

# Draw sizes mirror 03_rates_analysis/07_phylolm_jackknife.R so the whole-gene
# and subregion jackknives are directly comparable. The submitted subregion
# jackknife used a single 60% draw capped at 1000, which meant H3 (2187 usable
# branches) subsampled at 46% while H1 (612) subsampled at 60% — the two were
# not on the same footing. Each size now gets its own subdirectory.
DRAW_SIZES <- list(
  list(label = "n1000", fixed = 1000L, frac = NA),
  list(label = "pct15", fixed = NA,    frac = 0.15),
  list(label = "pct30", fixed = NA,    frac = 0.30),
  list(label = "pct60", fixed = NA,    frac = 0.60)
)
draws     <- if ("--draw" %in% args) args[which(args == "--draw") + 1] else
             vapply(DRAW_SIZES, `[[`, character(1), "label")
plot_draw <- if ("--plot-draw" %in% args) args[which(args == "--plot-draw") + 1] else "pct60"
force_run <- "--force" %in% args
min_n     <- if ("--min-n"    %in% args) as.integer(args[which(args == "--min-n")    + 1]) else 40L
base_seed <- if ("--seed"     %in% args) as.integer(args[which(args == "--seed")     + 1]) else 42L
n_thread  <- if ("--threads"  %in% args) as.integer(args[which(args == "--threads")  + 1]) else 4L
plot_only <- "--plot-only" %in% args

cap_blas_threads(n_thread)

PRED <- c("scaled_temp", "scaled_humidity",
          "scaled_delta_temp", "scaled_delta_humidity",
          "scaled_sampling_density", "scaled_pop_density")

RESPONSE_LABELS <- c(
  log_dN_rate = "dN (non-synonymous rate)",
  log_dS_rate = "dS (synonymous rate)",
  log_omega   = "omega (dN/dS)"
)

if (!coef_plot %in% PRED)
  warning("--coef '", coef_plot, "' is not one of: ", paste(PRED, collapse = ", "),
          immediate. = TRUE)

res_dir <- file.path(PATHS$rates, "jackknife", "slac")
fig_dir <- here("revisions", "figures", "jackknife", "slac")
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

stem_for <- function(gene, region, response)
  paste0(gene, "_", region, "_", response)

#' Directory for one draw size, e.g. .../jackknife/slac/pct60/
draw_dir <- function(label) file.path(res_dir, label)

#' Results path for one gene x region x response at one draw size.
res_path_for <- function(label, gene, region, response)
  file.path(draw_dir(label),
            paste0(stem_for(gene, region, response), "_results.csv"))

#' Deterministic per-combination seed, so re-running one region reproduces it.
combo_seed <- function(gene, region, response, label = "") {
  key <- paste(gene, region, response, label)
  base_seed + sum(utf8ToInt(key)) %% 100000L
}

#' Build the three SLAC responses for one region. Identical to 16_phylolm_slac.R.
add_responses <- function(d, r) {
  N  <- d[[paste0("slac_", r, "_N")]];  S  <- d[[paste0("slac_", r, "_S")]]
  EN <- d[[paste0("slac_", r, "_EN")]]; ES <- d[[paste0("slac_", r, "_ES")]]
  okN <- !is.na(N) & !is.na(EN) & EN > 0 & !is.na(d$t) & d$t > 0 & N > 0
  okS <- !is.na(S) & !is.na(ES) & ES > 0 & !is.na(d$t) & d$t > 0 & S > 0
  okO <- !is.na(N) & !is.na(S) & !is.na(EN) & !is.na(ES) &
         EN > 0 & ES > 0 & N > 0 & S > 0
  d$log_dN_rate <- ifelse(okN, log(N / (EN * d$t)), NA_real_)
  d$log_dS_rate <- ifelse(okS, log(S / (ES * d$t)), NA_real_)
  d$log_omega   <- ifelse(okO, log((N / EN) / (S / ES)), NA_real_)
  d
}

# ── Jackknife ─────────────────────────────────────────────────────────────────

if (!plot_only) for (gene in genes) {
  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_slac.csv"))
  tree_path  <- file.path(here(), "data", "trees", "timed",
                          paste0(gene, "_timetree.nwk"))
  if (!file.exists(rates_path)) {
    warning(gene, ": ", basename(rates_path), " not found — run 15_slac_dnds.R"); next
  }
  if (!file.exists(tree_path)) { warning(gene, ": timetree not found"); next }

  df        <- suppressWarnings(read_csv(rates_path, show_col_types = FALSE,
                                         guess_max = 50000))
  full_tree <- read.tree(tree_path)
  df$t      <- terminal_branch_times(full_tree)[df$accession]

  for (region in regions) {
    need <- paste0("slac_", region, c("_N", "_S", "_EN", "_ES"))
    if (!all(need %in% names(df))) {
      message("[", gene, " ", region, "] SLAC columns missing — skipping"); next
    }
    dr <- add_responses(df, region)

    for (response in responses) {
      tag <- paste0("[", gene, " ", region, " ", response, "]")
      if (!response %in% names(dr)) {
        message(tag, " unknown response — skipping"); next
      }

      keep <- intersect(c("accession", response, PRED, "year"), names(dr))
      sub  <- dr[, keep, drop = FALSE]
      sub  <- sub[stats::complete.cases(sub) & is.finite(sub[[response]]), ,
                  drop = FALSE]
      if (nrow(sub) < min_n) {
        message(tag, " only ", nrow(sub), " branches (< ", min_n, ") — skipping")
        next
      }

      prepped  <- prep_tree_data(sub, full_tree, quiet = TRUE)
      pool_df  <- prepped$data
      pool_tree <- prepped$tree

      n_avail <- nrow(pool_df)
      n_par   <- length(PRED) + 1L + length(unique(pool_df$year)) - 1L
      fml <- as.formula(paste(response, "~",
                              paste(c(PRED, "as.factor(year)"), collapse = " + ")))
      message("\n", tag, " ", n_avail, " branches available")

      for (ds in DRAW_SIZES) {
      if (!(ds$label %in% draws)) next

      out_path <- res_path_for(ds$label, gene, region, response)
      if (file.exists(out_path) && !force_run) {
        message("  [", ds$label, "] exists, skipping (use --force): ",
                basename(out_path))
        next
      }
      dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

      # Fixed sizes cap at what is available; fractions scale with it.
      size <- if (!is.na(ds$fixed)) min(ds$fixed, n_avail)
              else round(n_avail * ds$frac)

      # A draw the size of the whole pool is not a jackknife: sampling without
      # replacement just permutes every branch back in, so all replicates are
      # the same fit. This bites the fixed n1000 size on the smaller subregions
      # (H1 epitope has 612 usable branches, the ridge 125-187).
      if (size >= n_avail) {
        message("  [", ds$label, "] draw size ", size, " is the whole pool (",
                n_avail, ") — no subsampling possible, skipping")
        next
      }

      message("  [", ds$label, "] ", n_samples, " replicates of ", size,
              " (", round(100 * size / n_avail), "% of pool)  obs/par = ",
              round(size / n_par, 1))
      if (size / n_par < 10)
        message("    NOTE: under 10 observations per parameter — expect a wide ",
                "spread; that instability is itself the result")
      if (size < 20) { message("    replicate size < 20 — skipping"); next }

      set.seed(combo_seed(gene, region, response, ds$label))
      res_list <- vector("list", n_samples)
      mem_list <- vector("list", n_samples)
      n_fail   <- 0L
      t0 <- proc.time()

      for (i in seq_len(n_samples)) {
        acc  <- sample(pool_df$accession, size, replace = FALSE)
        sdf  <- pool_df[match(acc, pool_df$accession), , drop = FALSE]
        rownames(sdf) <- sdf$accession
        stree <- tryCatch(unroot(keep.tip(pool_tree, acc)),
                          error = function(e) NULL)
        if (is.null(stree)) { n_fail <- n_fail + 1L; next }

        fit <- tryCatch(
          phylolm(fml, data = sdf, phy = stree, model = PARAMS$phylolm_method),
          error = function(e) NULL)
        if (is.null(fit)) { n_fail <- n_fail + 1L; next }

        cm <- summary(fit)$coefficients
        se_col <- intersect(c("StdErr", "Std.Error", "std_error"),
                            colnames(cm))[1]
        t_col  <- intersect(c("t.value", "t value", "z.value"),
                            colnames(cm))[1]
        # Every predictor is stored, so --coef can change without refitting.
        terms_keep <- intersect(PRED, rownames(cm))
        res_list[[i]] <- tibble(
          gene = gene, region = region, response = response,
          sample_id = i, n = size,
          term      = terms_keep,
          estimate  = cm[terms_keep, "Estimate"],
          std_error = cm[terms_keep, se_col],
          t_value   = if (!is.na(t_col)) cm[terms_keep, t_col] else NA_real_,
          p.value   = cm[terms_keep, "p.value"],
          lambda    = tryCatch(fit$optpar, error = function(e) NA_real_)
        )
        mem_list[[i]] <- tibble(sample_id = i, accession = acc)

        if (i %% 200 == 0)
          message("  ", i, "/", n_samples, "  [",
                  round((proc.time() - t0)[["elapsed"]] / 60, 1), " min]")
      }

      results_df <- bind_rows(res_list)
      if (!nrow(results_df)) { message("  all replicates failed"); next }

      results_df$draw_size <- ds$label
      write_csv(results_df, out_path)
      write_csv(bind_rows(mem_list),
                sub("_results\\.csv$", "_membership.csv", out_path))

      ok <- results_df |> filter(term == coef_plot)
      message("    done in ", round((proc.time() - t0)[["elapsed"]] / 60, 1),
              " min; ", n_fail, " replicate(s) failed")
      if (nrow(ok))
        message("    ", coef_plot, ": median ", round(median(ok$estimate), 4),
                "  range [", round(min(ok$estimate), 4), ", ",
                round(max(ok$estimate), 4), "]  ",
                round(100 * mean(ok$p.value < 0.05), 1), "% of replicates p<0.05")
      }   # draw size
    }
  }
}

# ── Plots ─────────────────────────────────────────────────────────────────────

SIG_LEVELS  <- c("ns", "p < 0.05", "p < 0.01", "p < 0.001")
SIG_COLOURS <- setNames(viridis(4, end = 0.92), SIG_LEVELS)

classify_sig <- function(p) {
  factor(ifelse(is.na(p), "ns",
         ifelse(p < 0.001, "p < 0.001",
         ifelse(p < 0.01,  "p < 0.01",
         ifelse(p < 0.05,  "p < 0.05", "ns")))), levels = SIG_LEVELS)
}

summary_rows <- list()

for (gene in genes) for (region in regions) for (response in responses) {
  st    <- stem_for(gene, region, response)
  rpath <- res_path_for(plot_draw, gene, region, response)
  if (!file.exists(rpath)) next

  all_terms <- read_csv(rpath, show_col_types = FALSE)
  results   <- all_terms |> filter(term == coef_plot)
  if (!nrow(results)) {
    message("[", st, "] coefficient '", coef_plot, "' not in results — skipping")
    next
  }
  results$sig <- classify_sig(results$p.value)
  summary_rows[[st]] <- results

  sub_lab <- paste0(gene, "  |  ", REGION_LABELS[region], "  |  ",
                    RESPONSE_LABELS[[response]], "  |  n = ", results$n[1],
                    " per replicate")
  base <- file.path(fig_dir, paste0(st, "_", plot_draw, "_", coef_plot))

  # ── A: replicate order ──────────────────────────────────────────────────────
  pA <- ggplot(results, aes(x = sample_id, y = estimate, colour = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = estimate - std_error,
                       ymax = estimate + std_error),
                   linewidth = 0.3, alpha = 0.5) +
    geom_point(size = 1.2) +
    scale_colour_manual(values = SIG_COLOURS, name = NULL, drop = FALSE) +
    labs(title = paste0("Jackknife replicates: ", coef_plot),
         subtitle = sub_lab, x = "Replicate", y = "Coefficient") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(paste0(base, "_samples.png"), pA, width = 9, height = 4, dpi = 200)

  # ── B: ordered by estimate ──────────────────────────────────────────────────
  ordered <- results |> arrange(estimate) |> mutate(rank = row_number())
  pB <- ggplot(ordered, aes(x = rank, y = estimate, colour = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = estimate - std_error,
                       ymax = estimate + std_error),
                   linewidth = 0.3, alpha = 0.5) +
    geom_point(size = 1.2) +
    scale_colour_manual(values = SIG_COLOURS, name = NULL, drop = FALSE) +
    labs(title = paste0("Jackknife replicates, ordered: ", coef_plot),
         subtitle = sub_lab,
         x = "Replicate (ordered by estimate)", y = "Coefficient") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(paste0(base, "_samples_ordered.png"), pB, width = 9, height = 4, dpi = 200)

  # ── C: per-sequence ribbon ──────────────────────────────────────────────────
  # For each branch, the distribution of coefficients across the replicates it
  # appeared in. A sequence driving the result shows a visibly offset band.
  mpath <- sub("_results\\.csv$", "_membership.csv", rpath)
  if (file.exists(mpath)) {
    mem <- read_csv(mpath, show_col_types = FALSE)
    per_seq <- mem |>
      inner_join(results |> select(sample_id, estimate), by = "sample_id") |>
      group_by(accession) |>
      summarise(n_sampled = n(),
                median_est = median(estimate, na.rm = TRUE),
                q05 = quantile(estimate, 0.05, na.rm = TRUE),
                q25 = quantile(estimate, 0.25, na.rm = TRUE),
                q75 = quantile(estimate, 0.75, na.rm = TRUE),
                q95 = quantile(estimate, 0.95, na.rm = TRUE),
                .groups = "drop") |>
      arrange(median_est) |> mutate(rank = row_number())

    pC <- ggplot(per_seq, aes(x = rank)) +
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
      geom_ribbon(aes(ymin = q05, ymax = q95), alpha = 0.25,
                  fill = viridis(1, begin = 0.35)) +
      geom_ribbon(aes(ymin = q25, ymax = q75), alpha = 0.45,
                  fill = viridis(1, begin = 0.35)) +
      geom_line(aes(y = median_est), linewidth = 0.5) +
      labs(title = paste0("Per-sequence influence: ", coef_plot),
           subtitle = paste0(sub_lab, "  |  ribbons = IQR (dark) and 5-95% (light)"),
           x = "Sequence rank (ordered by median estimate)",
           y = "Coefficient across replicates containing that sequence") +
      theme_bw(base_size = 11) +
      theme(panel.grid.minor = element_blank())
    ggsave(paste0(base, "_sequences.png"), pC, width = 9, height = 4.5, dpi = 200)
  }
  message("  Saved figures: ", basename(base), "_*.png")
}

# ── Combined summary across regions and genes ─────────────────────────────────

if (length(summary_rows) > 0) {
  comb <- bind_rows(summary_rows) |>
    mutate(subtype   = ifelse(gene %in% c("H1", "N1"), "H1N1", "H3N2"),
           region_lab = factor(REGION_LABELS[region], levels = REGION_LEVELS)) |>
    group_by(gene, region, response) |>
    arrange(estimate, .by_group = TRUE) |>
    mutate(rank = row_number()) |>
    ungroup()

  annot <- comb |>
    group_by(subtype, region_lab) |>
    summarise(med = median(estimate, na.rm = TRUE),
              pct_sig = 100 * mean(p.value < 0.05, na.rm = TRUE),
              .groups = "drop") |>
    mutate(txt = sprintf("median %.3f\n%.0f%% p<0.05", med, pct_sig))

  p_sum <- ggplot(comb, aes(x = rank, y = estimate, colour = sig)) +
    geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
    geom_linerange(aes(ymin = estimate - std_error, ymax = estimate + std_error),
                   linewidth = 0.25, alpha = 0.4) +
    geom_point(size = 0.7, alpha = 0.8) +
    geom_text(data = annot, aes(x = Inf, y = Inf, label = txt),
              hjust = 1.05, vjust = 1.2, inherit.aes = FALSE,
              size = 2.6, colour = "grey25") +
    scale_colour_manual(values = SIG_COLOURS, name = NULL, drop = FALSE) +
    facet_grid(subtype ~ region_lab, scales = "free_y") +
    labs(title = paste0("Jackknife stability of ", coef_plot),
         subtitle = paste0(RESPONSE_LABELS[[responses[1]]],
                           "  |  replicates ordered by estimate within panel"),
         x = "Replicate (ordered by estimate)", y = "Coefficient") +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold", size = 9),
          panel.grid.minor = element_blank())

  out <- file.path(fig_dir, paste0("jackknife_summary_", responses[1], "_",
                                   plot_draw, "_", coef_plot, ".png"))
  ggsave(out, p_sum, width = 4 + 2.4 * n_distinct(comb$region), height = 6,
         dpi = 200, limitsize = FALSE)
  message("\nSaved: ", basename(out))

  message("\nSummary (", coef_plot, "):")
  comb |>
    group_by(gene, region, response) |>
    summarise(replicates = n(), n_per_rep = first(n),
              median = round(median(estimate), 4),
              q025 = round(quantile(estimate, 0.025), 4),
              q975 = round(quantile(estimate, 0.975), 4),
              pct_p05 = round(100 * mean(p.value < 0.05), 1),
              pct_same_sign = round(100 * mean(sign(estimate) ==
                                    sign(median(estimate))), 1),
              .groups = "drop") |>
    as.data.frame() |>
    print(row.names = FALSE)
}

message("\nDone.")
