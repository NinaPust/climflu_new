# scripts/03_rates_analysis/09_jackknife_sig_explore.R
#
# PURPOSE: Explore what distinguishes jackknife samples at a given significance
#          level. For each gene, filters samples to a chosen significance tier
#          and produces plots comparing significant vs non-significant samples:
#
#   1. Rate distribution
#   2. Geographic distribution (country frequency)
#   3. Year and month distributions
#   4. Temperature distribution
#   5. Sequence enrichment — sequences disproportionately present in
#      significant samples (potential signal drivers)
#
# INPUT:
#   data/rates/jackknife/{gene}/jackknife_results.csv
#   data/rates/jackknife/{gene}/jackknife_sequence_membership.csv
#   data/rates/{gene}_rates.csv
#
# OUTPUT:
#   figures/jackknife/explore/{gene}_{sig}/01_rates.png
#   figures/jackknife/explore/{gene}_{sig}/02_geography.png
#   figures/jackknife/explore/{gene}_{sig}/03_time.png
#   figures/jackknife/explore/{gene}_{sig}/04_temperature.png
#   figures/jackknife/explore/{gene}_{sig}/05_enrichment.png
#
# FLAGS:
#   --gene X   run a single gene (default: all GENES)
#   --sig X    significance tier to treat as "significant":
#              one of * / ** / *** (default: *)
#              Samples at this tier AND more significant are "significant";
#              ns-only samples are always the comparison baseline.

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(viridisLite)
  library(patchwork)
  library(maps)
})

source(here("config.R"))
source(here("R", "jackknife_utils.R"))
source(here("R", "climate_utils.R"))

# ── CLI args ──────────────────────────────────────────────────────────────────
args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) {
  args[which(args == "--gene") + 1]
} else {
  GENES
}
sig_arg <- if ("--sig" %in% args) {
  args[which(args == "--sig") + 1]
} else {
  "*"
}

SIG_ORDER <- c("ns", "*", "**", "***")
if (!sig_arg %in% SIG_ORDER) {
  stop("--sig must be one of: ns * ** *** (got: '", sig_arg, "')")
}

# Samples at sig_arg or more significant are treated as "significant"
sig_threshold_idx <- which(SIG_ORDER == sig_arg)
sig_tiers         <- SIG_ORDER[sig_threshold_idx:length(SIG_ORDER)]

classify_sig <- function(p) {
  factor(
    dplyr::case_when(
      p < 0.001 ~ "***",
      p < 0.01  ~ "**",
      p < 0.05  ~ "*",
      TRUE      ~ "ns"
    ),
    levels = SIG_ORDER
  )
}

TOP_N_ENRICHED <- 40L   # sequences to show in enrichment plot
JK_COEF <- "scaled_temp"

for (gene in genes) {
  message("\n[", gene, "] Exploring significance tier: ", sig_arg, " and above...")

  jk_dir     <- file.path(PATHS$rates, "jackknife", gene)
  res_path   <- file.path(jk_dir, "jackknife_results.csv")
  mem_path   <- file.path(jk_dir, "jackknife_sequence_membership.csv")
  rates_path <- file.path(paste0("/home/nina/climflu_new/data/rates/", gene, "_rates.csv"))

  for (p in c(res_path, mem_path, rates_path)) {
    if (!file.exists(p)) { warning("Missing: ", p); next }
  }

  results    <- read_jackknife(res_path, JK_COEF)
  membership <- read_csv(mem_path,    show_col_types = FALSE)
  rates      <- read_csv(rates_path,  show_col_types = FALSE)
  rates      <- scale_env_vars(rates)

  results <- results |>
    mutate(sig    = classify_sig(p.value),
           is_sig = sig %in% sig_tiers,
           is_ns  = sig == "ns",
           sig_label = case_when(
             sig %in% sig_tiers ~ paste0(sig_arg, "+"),
             sig == "ns"        ~ "ns",
             TRUE               ~ NA_character_
           ),
           sig_label = factor(
             sig_label, levels = c("ns", paste0(sig_arg, "+"))
           ))

  n_sig    <- sum(results$is_sig)
  n_ns     <- sum(results$is_ns)
  n_nonsig <- n_ns
  message(
    "  ", n_sig, " significant / ", n_ns,
    " ns samples (intermediate tiers excluded)"
  )

  if (n_sig == 0) {
    warning("  No samples at tier ", sig_arg, " or above — skipping ", gene)
    next
  }

  # ── Join sequences — keep only sig and ns samples (drop intermediate tiers)
  seq_sig <- membership |>
    left_join(results |> select(sample_id, is_sig, is_ns, sig_label, sig),
              by = "sample_id") |>
    filter(!is.na(sig_label)) |>
    left_join(rates, by = "accession")

  sig_label_pos <- paste0(sig_arg, "+")
  sig_colours   <- setNames(
    viridis(2, option = "D", begin = 0.2, end = 0.7),
    c("ns", sig_label_pos)
  )

  out_dir <- here("figures", "jackknife", "explore",
                  paste0(gene, "_", gsub("\\*", "star", sig_arg)))
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # ── 1. Rate distribution ───────────────────────────────────────────────────
  p1 <- ggplot(seq_sig |> filter(!is.na(log_rate)),
               aes(x = log_rate, fill = sig_label, colour = sig_label)) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    scale_fill_manual(
      values = sig_colours, drop = FALSE
    ) +
    scale_colour_manual(
      values = sig_colours, drop = FALSE
    ) +
    labs(
      title    = paste0(gene, " — rate distribution by significance"),
      x        = "log(rate)",
      y        = "Density",
      fill     = NULL, colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "01_rates.png"), p1,
         width = 7, height = 4, dpi = 180)

  # ── 2. Geographic distribution ────────────────────────────────────────────
  country_freq <- seq_sig |>
    filter(!is.na(country)) |>
    count(sig_label, country) |>
    group_by(sig_label) |>
    mutate(prop = n / sum(n)) |>
    ungroup()

  top_countries <- country_freq |>
    group_by(country) |>
    summarise(total = sum(n), .groups = "drop") |>
    slice_max(total, n = 20) |>
    pull(country)

  p2 <- country_freq |>
    filter(country %in% top_countries) |>
    mutate(country = reorder(country, prop)) |>
    ggplot(aes(x = prop, y = country, fill = sig_label)) +
    geom_col(position = "dodge", width = 0.7) +
    scale_fill_manual(
      values = sig_colours, drop = FALSE
    ) +
    scale_x_continuous(labels = scales::percent_format()) +
    labs(
      title = paste0(gene, " — geographic distribution (top 20 countries)"),
      x = "Proportion of sequences in sample group",
      y = NULL,
      fill = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "02_geography.png"), p2,
         width = 9, height = 7, dpi = 180)

  # ── 3. Year and month distributions ───────────────────────────────────────
  p3a <- seq_sig |>
    filter(!is.na(year)) |>
    count(sig_label, year) |>
    group_by(sig_label) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    ggplot(aes(x = year, y = prop, fill = sig_label, colour = sig_label)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
    scale_fill_manual(
      values = sig_colours, drop = FALSE
    ) +
    scale_colour_manual(
      values = sig_colours, drop = FALSE
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Year", y = "Proportion", fill = NULL, colour = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  p3b <- seq_sig |>
    filter(!is.na(month)) |>
    count(sig_label, month) |>
    group_by(sig_label) |>
    mutate(prop = n / sum(n)) |>
    ungroup() |>
    ggplot(aes(x = factor(month), y = prop,
               fill = sig_label, colour = sig_label)) +
    geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
    scale_fill_manual(
      values = sig_colours, drop = FALSE
    ) +
    scale_colour_manual(
      values = sig_colours, drop = FALSE
    ) +
    scale_x_discrete(
      labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                 "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
    ) +
    scale_y_continuous(labels = scales::percent_format()) +
    labs(x = "Month", y = "Proportion", fill = NULL, colour = NULL) +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  p3 <- (p3a / p3b) +
    plot_annotation(
      title = paste0(gene, " — temporal distribution by significance")
    )

  ggsave(file.path(out_dir, "03_time.png"), p3,
         width = 10, height = 7, dpi = 180)

  # ── 4. Temperature distribution ───────────────────────────────────────────
  p4 <- ggplot(seq_sig |> filter(!is.na(era5_temp_mean)),
               aes(x = era5_temp_mean, fill = sig_label,
                   colour = sig_label)) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    scale_fill_manual(
      values = sig_colours, drop = FALSE
    ) +
    scale_colour_manual(
      values = sig_colours, drop = FALSE
    ) +
    labs(
      title  = paste0(gene, " — temperature distribution by significance"),
      x      = "ERA5 temperature (°C)",
      y      = "Density",
      fill   = NULL, colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "04_temperature.png"), p4,
         width = 7, height = 4, dpi = 180)

  # ── 4b. Humidity distribution ─────────────────────────────────────────────
  p4b <- ggplot(seq_sig |> filter(!is.na(humidity_rh)),
                aes(x = humidity_rh, fill = sig_label, colour = sig_label)) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    scale_fill_manual(values = sig_colours, drop = FALSE) +
    scale_colour_manual(values = sig_colours, drop = FALSE) +
    labs(
      title  = paste0(gene, " — humidity distribution by significance"),
      x      = "Relative humidity (%)",
      y      = "Density",
      fill   = NULL, colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "04b_humidity.png"), p4b,
         width = 7, height = 4, dpi = 180)

  # ── 5. Latitude distribution ──────────────────────────────────────────────
  p5 <- ggplot(seq_sig |> filter(!is.na(latitude)),
               aes(x = latitude, fill = sig_label, colour = sig_label)) +
    geom_density(alpha = 0.35, linewidth = 0.5) +
    geom_vline(xintercept = 0, linetype = "dotted",
               colour = "grey40", linewidth = 0.4) +
    annotate("text", x = -2, y = Inf, label = "S", vjust = 1.5,
             hjust = 1, size = 3, colour = "grey40") +
    annotate("text", x = 2, y = Inf, label = "N", vjust = 1.5,
             hjust = 0, size = 3, colour = "grey40") +
    scale_fill_manual(values = sig_colours, drop = FALSE) +
    scale_colour_manual(values = sig_colours, drop = FALSE) +
    labs(
      title  = paste0(gene, " — latitude distribution by significance"),
      x      = "Latitude",
      y      = "Density",
      fill   = NULL, colour = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "05_latitude.png"), p5,
         width = 7, height = 4, dpi = 180)

  # ── 6. Sequence enrichment (raw ratio) ────────────────────────────────────
  seq_counts <- membership |>
    left_join(results |> select(sample_id, is_sig, is_ns),
              by = "sample_id") |>
    group_by(accession) |>
    summarise(
      n_sig_samples = sum(is_sig),
      n_ns_samples  = sum(is_ns),
      .groups       = "drop"
    ) |>
    mutate(
      freq_sig    = n_sig_samples / max(n_sig, 1),
      freq_nonsig = n_ns_samples  / max(n_ns,  1),
      ratio       = (freq_sig + 1e-6) / (freq_nonsig + 1e-6),
      n_total     = n_sig_samples + n_ns_samples
    ) |>
    arrange(desc(abs(ratio - 1)))

  top_enriched <- seq_counts |>
    slice_max(abs(ratio - 1), n = TOP_N_ENRICHED) |>
    mutate(accession = reorder(accession, ratio))

  p6 <- ggplot(top_enriched,
               aes(x = ratio, y = accession, fill = n_total)) +
    geom_vline(xintercept = 1, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_col(width = 0.7) +
    geom_text(
      aes(label = paste0(n_sig_samples, "/", n_sig,
                          " vs ", n_ns_samples, "/", n_ns)),
      hjust = -0.1, size = 2.2, colour = "grey30"
    ) +
    scale_fill_viridis_c(name = "Times\nsampled", option = "D") +
    labs(
      title    = paste0(
        gene, " — sequence enrichment in ", sig_arg,
        "+ vs ns samples (top ", TOP_N_ENRICHED, ")"
      ),
      subtitle = paste0(
        "freq in sig / freq in ns; ratio > 1 = more common in ",
        sig_arg, "+ samples\nthan in ns samples"
      ),
      x = paste0(
        "freq in ", sig_arg, "+ samples / freq in ns samples"
      ),
      y = NULL
    ) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    theme_bw(base_size = 9) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_blank(),
      axis.text.y        = element_text(size = 6)
    )

  ggsave(file.path(out_dir, "06_enrichment.png"), p6,
         width = 9, height = 8, dpi = 180)

  # ── 7. Sequence influence: median effect size vs times sampled ───────────
  seq_influence <- membership |>
    left_join(results |> select(sample_id, estimate, is_sig, is_ns),
              by = "sample_id") |>
    group_by(accession) |>
    summarise(
      n_total       = n(),
      n_sig_samples = sum(is_sig),
      n_ns_samples  = sum(is_ns),
      median_est    = median(estimate, na.rm = TRUE),
      mad_est       = mad(estimate, na.rm = TRUE),
      .groups       = "drop"
    ) |>
    mutate(prop_sig = n_sig_samples / n_total)

  p7 <- ggplot(
    seq_influence,
    aes(x = n_total, y = median_est, colour = prop_sig)
  ) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_colour_viridis_c(
      name   = paste0("Prop. samples\nthat were ", sig_arg, "+"),
      option = "D"
    ) +
    labs(
      title    = paste0(
        gene, " — per-sequence influence on temperature estimate"
      ),
      subtitle = paste0(
        "Each point is one sequence. Sequences sampled many times with\n",
        "high median estimate consistently co-occur with large effects."
      ),
      x = "Times sampled (out of 1000)",
      y = "Median temperature estimate across samples"
    ) +
    theme_bw(base_size = 11) +
    theme(panel.grid.minor = element_blank())

  ggsave(file.path(out_dir, "07_seq_influence.png"), p7,
         width = 8, height = 5, dpi = 180)

  write_csv(seq_counts, file.path(out_dir, "enrichment_scores.csv"))
  write_csv(seq_influence, file.path(out_dir, "seq_influence_scores.csv"))

  # ── 8–9. Enriched sequence deep-dive (ratio > 2) ──────────────────────────
  ENRICH_THRESHOLD <- 2

  enriched_acc <- seq_counts |>
    filter(ratio > ENRICH_THRESHOLD) |>
    pull(accession)

  n_enriched <- length(enriched_acc)
  message("  ", n_enriched, " sequences with enrichment ratio > ",
          ENRICH_THRESHOLD)

  if (n_enriched == 0) {
    message("  Skipping enriched-sequence plots (none above threshold)")
  } else {
    enriched_df <- rates |>
      filter(!is.na(latitude), !is.na(longitude)) |>
      left_join(
        seq_counts |> select(accession, ratio),
        by = "accession"
      ) |>
      mutate(enriched = accession %in% enriched_acc)

    # ── 8. World map ──────────────────────────────────────────────────────────
    world <- map_data("world")

    p8 <- ggplot() +
      geom_polygon(
        data = world,
        aes(x = long, y = lat, group = group),
        fill = "grey88", colour = "white", linewidth = 0.15
      ) +
      geom_point(
        data = enriched_df |> filter(!enriched),
        aes(x = longitude, y = latitude),
        colour = "grey60", size = 0.5, alpha = 0.3
      ) +
      geom_point(
        data = enriched_df |> filter(enriched),
        aes(x = longitude, y = latitude, colour = ratio),
        size = 1.8, alpha = 0.8
      ) +
      scale_colour_viridis_c(name = "Enrichment\nratio", option = "D") +
      coord_fixed(1.3, xlim = c(-180, 180), ylim = c(-60, 85)) +
      labs(
        title    = paste0(
          gene, " — location of enriched sequences (ratio > ",
          ENRICH_THRESHOLD, ")"
        ),
        subtitle = paste0(
          n_enriched, " enriched sequences (coloured) vs all ",
          "sequences in dataset (grey)"
        ),
        x = NULL, y = NULL
      ) +
      theme_bw(base_size = 11) +
      theme(
        panel.grid       = element_blank(),
        axis.text        = element_blank(),
        axis.ticks       = element_blank()
      )

    ggsave(file.path(out_dir, "08_enriched_map.png"), p8,
           width = 12, height = 6, dpi = 180)

    # ── 9. Enriched sequence distributions (rate/time/temp/humidity) ──────────
    enrich_label  <- paste0("enriched (ratio>", ENRICH_THRESHOLD, ")")
    enriched_long <- rates |>
      mutate(group = ifelse(
        accession %in% enriched_acc, enrich_label, "background"
      ),
      group = factor(group, levels = c("background", enrich_label)))

    dist_colours <- setNames(
      viridis(2, option = "D", begin = 0.2, end = 0.7),
      c("background", enrich_label)
    )

    make_density <- function(data, var, xlab) {
      ggplot(data |> filter(!is.na(.data[[var]])),
             aes(x = .data[[var]], fill = group, colour = group)) +
        geom_density(alpha = 0.35, linewidth = 0.5) +
        scale_fill_manual(values = dist_colours, drop = FALSE) +
        scale_colour_manual(values = dist_colours, drop = FALSE) +
        labs(x = xlab, y = "Density", fill = NULL, colour = NULL) +
        theme_bw(base_size = 10) +
        theme(panel.grid.minor = element_blank())
    }

    p9a <- make_density(enriched_long, "log_rate", "log(rate)")
    p9b <- make_density(enriched_long, "era5_temp_mean", "Temperature (°C)")
    p9c <- make_density(enriched_long, "humidity_rh", "Relative humidity (%)")

    p9d <- enriched_long |>
      filter(!is.na(year)) |>
      count(group, year) |>
      group_by(group) |>
      mutate(prop = n / sum(n)) |>
      ungroup() |>
      ggplot(aes(x = year, y = prop, fill = group, colour = group)) +
      geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
      scale_fill_manual(values = dist_colours, drop = FALSE) +
      scale_colour_manual(values = dist_colours, drop = FALSE) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "Year", y = "Proportion", fill = NULL, colour = NULL) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())

    p9e <- enriched_long |>
      filter(!is.na(month)) |>
      count(group, month) |>
      group_by(group) |>
      mutate(prop = n / sum(n)) |>
      ungroup() |>
      ggplot(aes(x = factor(month), y = prop,
                 fill = group, colour = group)) +
      geom_col(position = "dodge", width = 0.7, alpha = 0.8) +
      scale_fill_manual(values = dist_colours, drop = FALSE) +
      scale_colour_manual(values = dist_colours, drop = FALSE) +
      scale_x_discrete(
        labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun",
                   "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
      ) +
      scale_y_continuous(labels = scales::percent_format()) +
      labs(x = "Month", y = "Proportion", fill = NULL, colour = NULL) +
      theme_bw(base_size = 10) +
      theme(panel.grid.minor = element_blank())

    p9 <- (p9a | p9b | p9c) / (p9d | p9e) +
      plot_annotation(
        title    = paste0(
          gene, " — enriched sequence distributions (ratio > ",
          ENRICH_THRESHOLD, ")"
        ),
        subtitle = paste0(
          n_enriched, " enriched vs ",
          nrow(rates) - n_enriched, " background sequences"
        )
      ) +
      plot_layout(guides = "collect")

    ggsave(file.path(out_dir, "09_enriched_distributions.png"), p9,
           width = 14, height = 8, dpi = 180)
  }

  message("  Written: ", out_dir)
}

message("\nDone.")
