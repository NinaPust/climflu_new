# scripts/03_rates_analysis/05_inla_diagnostics.R
#
# PURPOSE: Diagnostic figures for the INLA subsample used in 04_inla_phylo.R.
#          For each gene, produces a multi-page PDF with:
#            Page 1 — Subsample tree (fan layout, coloured by year)
#            Page 2 — Density curves for log_rate
#            Page 3 — Density curves for climate/population predictors
#          Full dataset shown as grey fill; subsample as coloured line.
#          KS tests (subsample vs full) are printed to console and captioned
#          on each plot so you can spot distributional mismatch quantitatively.
#
# INPUT:
#   data/rates/{gene}_rates.csv
#   data/trees/timed/{gene}_timetree.nwk
#
# OUTPUT:
#   figures/inla_diagnostics/{gene}_2000.pdf
#
# FLAGS:
#   --gene X   run a single gene segment

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

source(here("config.R"))

ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")
args      <- commandArgs(trailingOnly = TRUE)
genes     <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  ALL_GENES
}

N_SUB <- 2000

out_dir <- here("figures", "inla_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Helpers ───────────────────────────────────────────────────────────────────

stratified_subsample <- function(df, n, seed = 42) {
  set.seed(seed)
  if (nrow(df) <= n) return(df)
  total_n     <- nrow(df)
  year_counts <- df |>
    count(.data$year, name = "grp_n") |>
    mutate(.k = pmax(1L, as.integer(round(.data$grp_n / total_n * n))))
  sampled <- lapply(split(df, df$year), function(sub) {
    k <- year_counts$.k[year_counts$year == sub$year[1]]
    sub[sample(nrow(sub), min(k, nrow(sub))), ]
  })
  result <- do.call(rbind, sampled)
  result[sample(nrow(result), min(n, nrow(result))), ]
}

# Density overlay: full dataset as grey fill, subsample as coloured line.
# Returns list(plot, ks) so the caller can assemble captions.
overlay_density <- function(full_df, sub_df, col, label, line_col = "#2166ac") {
  vals_full <- full_df[[col]]
  vals_sub  <- sub_df[[col]]
  vals_full <- vals_full[is.finite(vals_full)]
  vals_sub  <- vals_sub[is.finite(vals_sub)]
  if (length(vals_full) < 2 || length(vals_sub) < 2) return(NULL)

  ks  <- ks.test(vals_sub, vals_full)
  sig <- dplyr::case_when(
    ks$p.value < 0.001 ~ "***",
    ks$p.value < 0.01  ~ "**",
    ks$p.value < 0.05  ~ "*",
    TRUE               ~ "ns"
  )
  cap <- sprintf("KS D = %.3f, p = %.3f %s", ks$statistic, ks$p.value, sig)

  p <- ggplot() +
    geom_density(
      data  = data.frame(x = vals_full),
      aes(x = x),
      fill  = "grey85", color = "grey60", linewidth = 0.4
    ) +
    geom_density(
      data      = data.frame(x = vals_sub),
      aes(x = x),
      color     = line_col, fill = NA, linewidth = 0.9
    ) +
    labs(x = label, y = "Density", caption = cap) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor = element_blank(),
      plot.caption     = element_text(size = 7, color = "grey40")
    )

  list(plot = p, ks = ks, sig = sig, label = label)
}

# ── Main loop ─────────────────────────────────────────────────────────────────

for (gene in genes) {
  message("\n[", gene, "] Building diagnostics...")

  rates_path <- file.path(PATHS$rates, paste0(gene, "_rates.csv"))
  tree_path  <- file.path(PATHS$trees_timed, paste0(gene, "_timetree.nwk"))

  if (!file.exists(rates_path)) { warning("Missing: ", rates_path); next }
  if (!file.exists(tree_path))  { warning("Missing: ", tree_path);  next }

  full_df <- read_csv(rates_path, show_col_types = FALSE)
  tree    <- read.tree(tree_path)

  # Mirror 04_inla_phylo.R cleaning exactly
  common <- intersect(tree$tip.label, full_df$accession)
  tree   <- keep.tip(tree, common)

  tip_edge_idx <- match(seq_len(Ntip(tree)), tree$edge[, 2])
  zero_tips    <- tree$tip.label[tree$edge.length[tip_edge_idx] == 0]
  if (length(zero_tips) > 0) tree <- drop.tip(tree, zero_tips)
  common <- tree$tip.label

  full_df <- full_df |>
    filter(.data$accession %in% common, is.finite(.data$log_rate)) |>
    arrange(match(.data$accession, tree$tip.label))

  if ("scaled_pop_density" %in% names(full_df))
    full_df <- full_df |> filter(!is.na(.data$scaled_pop_density))

  tree <- keep.tip(tree, full_df$accession)

  sub_df   <- stratified_subsample(full_df, N_SUB)
  sub_tree <- keep.tip(tree, sub_df$accession)

  n_full <- nrow(full_df)
  n_sub  <- nrow(sub_df)
  message("  Full: ", n_full, "  Subsample: ", n_sub)

  # ── KS table printed to console ───────────────────────────────────────────
  all_vars <- c("log_rate", "scaled_temp", "scaled_humidity",
                "scaled_pop_density", "scaled_season")
  ks_rows <- lapply(all_vars, function(v) {
    if (!v %in% names(full_df)) return(NULL)
    vf <- full_df[[v]][is.finite(full_df[[v]])]
    vs <- sub_df[[v]][is.finite(sub_df[[v]])]
    if (length(vf) < 2 || length(vs) < 2) return(NULL)
    ks <- ks.test(vs, vf)
    data.frame(variable = v,
               D        = round(ks$statistic, 4),
               p        = round(ks$p.value,   4),
               stringsAsFactors = FALSE)
  })
  ks_table <- do.call(rbind, ks_rows)
  message("  KS tests (subsample vs full):")
  print(ks_table, row.names = FALSE)

  out_path <- file.path(out_dir, paste0(gene, "_2000.pdf"))
  pdf(out_path, width = 8.27, height = 11.69)

  # ── Page 1: Tree coloured by year ─────────────────────────────────────────
  tip_meta <- data.frame(
    label = sub_df$accession,
    year  = sub_df$year,
    stringsAsFactors = FALSE
  )

  p_tree <- ggtree(sub_tree, layout = "fan", open.angle = 15,
                   linewidth = 0.15, color = "grey60") %<+% tip_meta +
    geom_tippoint(aes(color = year), size = 0.4, shape = 16) +
    scale_color_viridis_c(option = "plasma", name = "Year") +
    labs(
      title    = paste0(gene, " — INLA subsample tree  (n = ", n_sub, ")"),
      subtitle = paste0("Full dataset after cleaning: n = ", n_full,
                        "  |  Subsample seed = 42, stratified by year")
    ) +
    theme_tree() +
    theme(
      plot.title      = element_text(size = 11, face = "bold"),
      plot.subtitle   = element_text(size = 8,  color = "grey40"),
      legend.position    = "right",
      legend.key.height  = unit(0.8, "cm")
    )

  print(p_tree)

  # ── Page 2: log_rate density ───────────────────────────────────────────────
  res_rate <- overlay_density(full_df, sub_df, "log_rate",
                              "log(substitution rate)", line_col = "#d73027")

  if (!is.null(res_rate)) {
    p_rate <- res_rate$plot +
      labs(
        title    = paste0(gene, " — log_rate distribution"),
        subtitle = "Grey fill = full dataset  |  Red line = INLA subsample"
      ) +
      theme(
        plot.title    = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8,  color = "grey40")
      )
    print(p_rate)
  }

  # ── Page 3: Climate + population predictor densities ──────────────────────
  predictor_meta <- list(
    scaled_temp        = list(label = "Scaled temperature",        col = "#e08214"),
    scaled_humidity    = list(label = "Scaled relative humidity",  col = "#542788"),
    scaled_pop_density = list(label = "Scaled population density", col = "#1b7837"),
    scaled_season      = list(label = "Scaled season (cosine)",    col = "#2166ac")
  )

  plots <- list()
  for (pred in names(predictor_meta)) {
    if (!pred %in% names(full_df)) next
    pm  <- predictor_meta[[pred]]
    res <- overlay_density(full_df, sub_df, pred, pm$label, line_col = pm$col)
    if (!is.null(res)) plots[[pred]] <- res$plot
  }

  if (length(plots) > 0) {
    combo <- Reduce(`+`, plots) +
      plot_layout(ncol = 2L) +
      plot_annotation(
        title    = paste0(gene, " — predictor distributions"),
        subtitle = "Grey fill = full dataset  |  Coloured line = INLA subsample  |  KS p-value captioned per panel",
        theme    = theme(
          plot.title    = element_text(size = 11, face = "bold"),
          plot.subtitle = element_text(size = 8,  color = "grey40")
        )
      )
    print(combo)
  }

  dev.off()
  message("  Saved: ", out_path)
}

message("\nDone.")
