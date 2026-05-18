# scripts/04_epitope_analysis/12_phylolm_ns_rates.R
#
# PURPOSE: Phylogenetic regression of NS/S substitution rates (from script 09)
#          against climate predictors for all three regions (ridge, epitope,
#          full gene) simultaneously. Avoids the dS denominator problem of
#          yn00-based omega by using raw substitution counts / branch length.
#
#          After fitting, produces:
#            figures/epitope/phylolm_forest_ns.png      (all predictors, p<0.05)
#            figures/epitope/phylolm_forest_ns_fdr.png  (BH-FDR p<0.05)
#            figures/epitope/{region}/ns_rate_density.png
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   data/rates/{gene}_phylolm_ns_results.csv
#   data/rates/all_phylolm_ns_results.csv
#   figures/epitope/phylolm_forest_ns.png
#   figures/epitope/phylolm_forest_ns_fdr.png
#   figures/epitope/{region}/ns_rate_density.png
#
# FLAGS:
#   --gene X   run a single gene (default: H1 H3)

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

source(here("config.R"))
source(here("R", "phylolm_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else
         c("H1", "H3")

predictor_terms <- c(
  "scaled_temp", "scaled_humidity",
  "scaled_delta_temp", "scaled_delta_humidity",
  "scaled_sampling_density", "scaled_pop_density"
)

responses <- c(
  "log_ridge_ns_rate",   "log_ridge_s_rate",
  "log_epitope_ns_rate", "log_epitope_s_rate",
  "log_full_ns_rate",    "log_full_s_rate",
  "log_HA2_ns_rate",     "log_HA2_s_rate"
)

predictor_labels <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative humidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔHumidity",
  scaled_sampling_density = "Sampling density",
  scaled_pop_density      = "Population density"
)

response_labels <- c(
  log_ridge_ns_rate   = "Ridge NS rate",
  log_ridge_s_rate    = "Ridge S rate",
  log_epitope_ns_rate = "Epitope NS rate",
  log_epitope_s_rate  = "Epitope S rate",
  log_full_ns_rate    = "Full gene NS rate",
  log_full_s_rate     = "Full gene S rate",
  log_HA2_ns_rate     = "HA2 NS rate",
  log_HA2_s_rate      = "HA2 S rate"
)

region_labels <- c(
  ridge   = "Ridge",
  epitope = "Epitope",
  full    = "Full gene",
  HA2     = "HA2"
)

# ── Phylolm per gene ──────────────────────────────────────────────────────────

for (gene in genes) {
  message("\n[", gene, "] phylolm NS-rate analysis...")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates_epitope.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))
  out_path   <- file.path(PATHS$rates_phylolm, paste0(gene, "_phylolm_ns_results.csv"))

  if (!file.exists(rates_path)) { warning("Rates not found: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Tree not found: ",  tree_path);  next }

  df   <- read_csv(rates_path, show_col_types = FALSE)
  tree <- read.tree(tree_path)

  avail <- intersect(responses, names(df))
  avail <- avail[sapply(avail, function(r) sum(!is.na(df[[r]])) > 10)]
  if (length(avail) == 0) { message("  No response columns, skipping."); next }

  common <- intersect(tree$tip.label, df$accession)
  tree   <- keep.tip(tree, common)

  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) {
    message("  Dropping ", length(zero_tips), " zero-length tips")
    tree   <- drop.tip(tree, zero_tips)
    common <- tree$tip.label
  }

  df <- df |> filter(accession %in% common) |>
    arrange(match(accession, tree$tip.label))
  rownames(df) <- df$accession

  predictors <- intersect(predictor_terms, names(df))
  message("  ", nrow(df), " taxa  |  responses: ", paste(avail, collapse = ", "))

  all_results <- list()
  for (response in avail) {
    sub_df   <- df |> filter(!is.na(.data[[response]]))
    sub_tree <- keep.tip(tree, sub_df$accession)
    if (nrow(sub_df) < 10) next

    year_term    <- if ("year" %in% names(sub_df)) "as.factor(year)" else NULL
    formula_full <- as.formula(
      paste(response, "~",
            paste(c(predictors, year_term), collapse = " + "))
    )
    result <- run_phylolm(
      formula  = formula_full,
      data     = sub_df,
      tree     = sub_tree,
      method   = PARAMS$phylolm_method,
      subset   = "full",
      response = response,
      gene     = gene
    )
    all_results[[paste(gene, response, sep = "_")]] <- result
  }

  results_df <- bind_rows(all_results)
  write_csv(results_df, out_path)
  message("  Written: ", out_path)
  print(results_df |>
          filter(term %in% predictor_terms) |>
          select(gene, response, term, estimate, p.value))
}

# ── BH-FDR across all genes ───────────────────────────────────────────────────
message("\nApplying BH-FDR correction...")

all_res <- lapply(genes, function(g) {
  path <- file.path(PATHS$rates_phylolm, paste0(g, "_phylolm_ns_results.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE)
}) |> bind_rows()

if (nrow(all_res) == 0) stop("No phylolm results found.")

pred_rows     <- all_res$term %in% predictor_terms
all_res$p.fdr <- NA_real_
all_res$p.fdr[pred_rows] <- p.adjust(all_res$p.value[pred_rows], method = "BH")

combined_path <- file.path(PATHS$rates_phylolm, "all_phylolm_ns_results.csv")
write_csv(all_res, combined_path)
message("Written: ", combined_path)

sig <- all_res |>
  filter(term %in% predictor_terms, p.fdr < 0.05) |>
  select(gene, response, term, estimate, p.value, p.fdr) |>
  mutate(across(estimate, \(x) round(x, 4)),
         across(c(p.value, p.fdr), \(x) signif(x, 3))) |>
  arrange(p.fdr)
message("\nSignificant after FDR (p.fdr < 0.05):")
print(sig)

# ── Plots ─────────────────────────────────────────────────────────────────────

out_root <- here("figures", "epitope")
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

res <- all_res

pred_colours <- setNames(
  viridis(length(predictor_labels), end = 0.92),
  predictor_labels
)

plot_df <- res |>
  filter(term %in% names(predictor_labels),
         response %in% names(response_labels)) |>
  mutate(
    ci_lo     = estimate - 1.96 * std_error,
    ci_hi     = estimate + 1.96 * std_error,
    predictor = factor(predictor_labels[term], levels = predictor_labels),
    region    = factor(
      case_when(
        grepl("^log_ridge",   response) ~ "Ridge",
        grepl("^log_epitope", response) ~ "Epitope",
        grepl("^log_full",    response) ~ "Full gene",
        grepl("^log_HA2",     response) ~ "HA2"
      ),
      levels = c("Ridge", "Epitope", "Full gene", "HA2")
    ),
    rate_type = factor(
      ifelse(grepl("_ns_rate$", response), "NS", "S"),
      levels = c("NS", "S"),
      labels = c("NS rate", "S rate")
    ),
    subtype = ifelse(gene %in% c("H1", "N1"), "H1N1", "H3N2")
  )

make_forest <- function(df, sig_col, subtitle, filename) {
  df <- df |>
    mutate(point_fill = ifelse(.data[[sig_col]],
                               pred_colours[as.character(predictor)], "white"))
  p <- ggplot(df, aes(x = estimate, y = region,
                      colour = predictor, group = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.7,
                   position = position_dodge(width = 0.7)) +
    geom_point(aes(fill = point_fill), shape = 21, size = 3, stroke = 0.8,
               position = position_dodge(width = 0.7)) +
    scale_colour_manual(values = pred_colours, name = NULL) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev) +
    facet_grid(rate_type ~ subtype) +
    labs(
      title    = "Climate predictors of NS and S substitution rates",
      subtitle = subtitle,
      x = "Coefficient estimate", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(strip.text         = element_text(face = "bold", size = 10),
          legend.position    = "bottom",
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank())
  ggsave(file.path(out_root, filename), p, width = 10, height = 8, dpi = 180)
  message("Saved: ", filename)
}

make_forest(
  plot_df |> mutate(sig = !is.na(p.value) & p.value < 0.05),
  "sig",
  "Coefficients ± 95% CI  |  Filled = p < 0.05",
  "phylolm_forest_ns.png"
)
make_forest(
  plot_df |> mutate(sig = !is.na(p.fdr) & p.fdr < 0.05),
  "sig",
  "Coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05",
  "phylolm_forest_ns_fdr.png"
)

# ── NS rate density by region, saved per region subdir ───────────────────────
rate_cols <- list(
  ridge   = c(NS = "log_ridge_ns_rate",   S = "log_ridge_s_rate"),
  epitope = c(NS = "log_epitope_ns_rate", S = "log_epitope_s_rate"),
  full    = c(NS = "log_full_ns_rate",    S = "log_full_s_rate"),
  HA2     = c(NS = "log_HA2_ns_rate",     S = "log_HA2_s_rate")
)

for (reg in names(rate_cols)) {
  reg_dir <- here("figures", "epitope", if (reg == "full") "whole" else reg)
  dir.create(reg_dir, recursive = TRUE, showWarnings = FALSE)

  dat <- lapply(c("H1", "H3"), function(g) {
    path <- file.path(PATHS$rates, paste0(g, "_rates_epitope.csv"))
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE) |>
      mutate(gene = g, subtype = ifelse(g == "H1", "H1N1", "H3N2"))
  }) |> bind_rows()

  ns_col <- rate_cols[[reg]][["NS"]]
  s_col  <- rate_cols[[reg]][["S"]]
  if (!ns_col %in% names(dat)) next

  long <- dat |>
    select(gene, subtype,
           NS = !!ns_col, S = !!s_col) |>
    pivot_longer(c(NS, S), names_to = "rate_type", values_to = "log_rate") |>
    filter(is.finite(log_rate)) |>
    mutate(rate_type = factor(rate_type, levels = c("NS", "S")))

  reg_colours <- c(NS = viridis(3, option = "turbo")[1],
                   S  = viridis(3, option = "turbo")[3])

  p <- ggplot(long, aes(x = log_rate, fill = rate_type, colour = rate_type)) +
    geom_density(alpha = 0.4, linewidth = 0.7) +
    facet_wrap(~subtype) +
    scale_fill_manual(values = reg_colours, name = NULL) +
    scale_colour_manual(values = reg_colours, name = NULL) +
    labs(
      title = paste0("Distribution of ", tools::toTitleCase(reg), " NS/S rates"),
      x = "log(rate)", y = "Density"
    ) +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  ggsave(file.path(reg_dir, "ns_rate_density.png"),
         p, width = 9, height = 5, dpi = 180)
  message("Saved: ", reg, "/ns_rate_density.png")
}

message("\nDone.")
