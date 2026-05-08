# scripts/04_epitope_analysis/10_plot_substitutions.R
#
# PURPOSE: Visualise NS/S substitution counts and individual mutation events
#          from script 09 for any region (ridge / epitope / full).
#
#          Per-region figures are saved to figures/epitope/{region}/:
#            ns_counts.png        — NS count distribution per tip
#            ns_temporal.png      — NS fraction by year
#            ns_climate.png       — climate exposure NS=0 vs NS>0
#            mut_changes.png      — convergence bar (ridge/epitope only)
#            mut_temporal.png     — temporal stacked bar (ridge/epitope only)
#            mut_chord.png        — AA→AA chord diagram (ridge/epitope only)
#            mut_chord_nt.png     — nt→nt chord diagram (ridge/epitope only)
#            mut_positions.png    — mutation count by position (ridge only)
#
# INPUT:
#   data/rates/{gene}_rates_epitope.csv
#   data/epitope/{gene}_{region}_mutations.csv   (ridge/epitope)
#
# FLAGS:
#   --region X   ridge | epitope | full  (default: ridge epitope full)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(forcats)
  library(viridisLite)
  library(circlize)
})

source(here("config.R"))

args    <- commandArgs(trailingOnly = TRUE)
regions <- if ("--region" %in% args) args[which(args == "--region") + 1] else
           c("ridge", "epitope", "full")

# ── Helpers ───────────────────────────────────────────────────────────────────

load_rates <- function(genes = c("H1", "H3")) {
  lapply(genes, function(g) {
    path <- file.path(PATHS$rates, paste0(g, "_rates_epitope.csv"))
    if (!file.exists(path)) return(NULL)
    read_csv(path, show_col_types = FALSE) |>
      mutate(gene = g, subtype = ifelse(g == "H1", "H1N1", "H3N2"))
  }) |> bind_rows()
}

load_mutations <- function(region, genes = c("H1", "H3")) {
  lapply(genes, function(g) {
    path <- file.path(here("data", "epitope"),
                      paste0(g, "_", region, "_mutations.csv"))
    if (!file.exists(path)) return(NULL)
    rates_meta <- read_csv(
      file.path(PATHS$rates, paste0(g, "_rates_epitope.csv")),
      show_col_types = FALSE
    ) |> select(accession, year)
    read_csv(path, show_col_types = FALSE) |>
      left_join(rates_meta, by = "accession") |>
      mutate(subtype = ifelse(g == "H1", "H1N1", "H3N2"))
  }) |> bind_rows()
}

# ── Main loop ─────────────────────────────────────────────────────────────────

rates_all <- load_rates()

for (region in regions) {
  message("\n[", region, "] Plotting substitutions...")

  out_dir <- here("figures", "epitope", region)
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  pfx  <- region
  ns_col  <- paste0(pfx, "_ns_count")
  s_col   <- paste0(pfx, "_s_count")

  if (!ns_col %in% names(rates_all)) {
    message("  Column ", ns_col, " not found — run script 09 first. Skipping.")
    next
  }

  dat <- rates_all |>
    rename(ns_count = !!ns_col, s_count = !!s_col) |>
    filter(!is.na(ns_count))

  pal3 <- viridis(3, end = 0.85)

  # ── 1. NS count distribution ────────────────────────────────────────────────
  count_df <- dat |>
    count(subtype, ns_count) |>
    mutate(ns_count = factor(ns_count))

  p1 <- ggplot(count_df, aes(x = ns_count, y = n, fill = ns_count)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = n), vjust = -0.4, size = 3.5) +
    facet_wrap(~subtype) +
    scale_fill_viridis_d(end = 0.85) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title    = paste0(tools::toTitleCase(region), " NS substitution counts per tip"),
      subtitle = "Counts from ancestral parsimony reconstruction",
      x = "NS substitutions on terminal branch", y = "Number of tips"
    ) +
    theme_bw(base_size = 12)

  ggsave(file.path(out_dir, "ns_counts.png"),
         p1, width = 8, height = 5, dpi = 180)
  message("  Saved: ns_counts.png")

  # ── 2. Temporal distribution ─────────────────────────────────────────────────
  temporal_df <- dat |>
    filter(!is.na(year)) |>
    group_by(subtype, year) |>
    summarise(total_tips = n(), ns_tips = sum(ns_count > 0),
              ns_frac = ns_tips / total_tips, .groups = "drop")

  p2 <- ggplot(temporal_df, aes(x = year)) +
    geom_col(aes(y = total_tips), fill = "grey80", width = 0.8) +
    geom_col(aes(y = ns_tips), fill = pal3[3], width = 0.8) +
    geom_line(aes(y = ns_frac * max(total_tips)),
              colour = pal3[2], linewidth = 1) +
    geom_point(aes(y = ns_frac * max(total_tips)),
               colour = pal3[2], size = 2) +
    facet_wrap(~subtype, scales = "free_y") +
    scale_y_continuous(
      name = "Tips (bars)",
      sec.axis = sec_axis(~ . / max(temporal_df$total_tips),
                          name = "Fraction with NS > 0 (line)")
    ) +
    labs(
      title    = paste0(tools::toTitleCase(region), " NS substitutions by year"),
      subtitle = "Grey = all tips; yellow = tips with NS > 0; line = fraction",
      x = "Year"
    ) +
    theme_bw(base_size = 12) +
    theme(axis.title.y.right = element_text(colour = pal3[2]))

  ggsave(file.path(out_dir, "ns_temporal.png"),
         p2, width = 10, height = 5, dpi = 180)
  message("  Saved: ns_temporal.png")

  # ── 3. Climate exposure NS=0 vs NS>0 ─────────────────────────────────────────
  clim_vars <- c("era5_temp_mean", "humidity_rh", "delta_temp", "delta_humidity")
  clim_vars <- intersect(clim_vars, names(dat))

  if (length(clim_vars) > 0) {
    clim_labels <- c(era5_temp_mean = "Temperature", humidity_rh = "Relative humidity",
                     delta_temp = "ΔTemperature", delta_humidity = "ΔHumidity")
    climate_df <- dat |>
      filter(if_all(all_of(clim_vars), is.finite)) |>
      mutate(has_ns = factor(ifelse(ns_count > 0, "NS > 0", "NS = 0"),
                             levels = c("NS = 0", "NS > 0"))) |>
      select(subtype, has_ns, all_of(clim_vars)) |>
      rename(any_of(clim_labels)) |>
      pivot_longer(-c(subtype, has_ns),
                   names_to = "variable", values_to = "value")

    p3 <- ggplot(climate_df, aes(x = has_ns, y = value, fill = has_ns)) +
      geom_violin(alpha = 0.6, linewidth = 0.4, draw_quantiles = 0.5) +
      geom_boxplot(width = 0.12, outlier.shape = NA,
                   fill = "white", linewidth = 0.5) +
      facet_grid(variable ~ subtype, scales = "free_y") +
      scale_fill_manual(values = c("NS = 0" = pal3[1], "NS > 0" = pal3[3]),
                        name = NULL) +
      labs(
        title    = paste0("Climate exposure: tips with vs without ",
                          region, " NS substitution"),
        subtitle = "Violin + median boxplot",
        x = NULL, y = "Value"
      ) +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom", strip.text = element_text(size = 9))

    ggsave(file.path(out_dir, "ns_climate.png"),
           p3, width = 9, height = 9, dpi = 180)
    message("  Saved: ns_climate.png")
  }

  # ── Mutation-level plots (ridge and epitope only) ────────────────────────────
  if (!region %in% c("ridge", "epitope")) next

  mut_all <- load_mutations(region)
  if (is.null(mut_all) || nrow(mut_all) == 0) {
    message("  No mutations data — skipping mutation plots.")
    next
  }

  # For epitope: join site_map to get antigenic region label
  if (region == "epitope") {
    site_maps <- lapply(c("H1", "H3"), function(g) {
      sm_path <- file.path(here("data", "epitope"), paste0(g, "_site_map.csv"))
      if (!file.exists(sm_path)) return(NULL)
      read_csv(sm_path, show_col_types = FALSE) |>
        filter(grepl("^epitope", region)) |>
        mutate(epi_region = sub("_ridge$", "", region)) |>
        distinct(residue_pos, epi_region) |>
        mutate(gene = g)
    }) |> bind_rows()
    mut_all <- mut_all |>
      left_join(site_maps |> select(gene, residue_pos, epi_region),
                by = c("gene", "residue_pos"))
  }

  min_count <- 5L
  recurrent <- mut_all |>
    count(subtype, change) |>
    filter(n > min_count) |>
    select(subtype, change)

  mut_rec <- mut_all |> semi_join(recurrent, by = c("subtype", "change"))

  # ── 4. Ridge-only: mutation count by position ──────────────────────────────
  if (region == "ridge") {
    pos_df <- mut_all |>
      count(subtype, residue_pos) |>
      group_by(subtype) |>
      mutate(pct = n / sum(n) * 100,
             residue_pos = factor(residue_pos))

    p4 <- ggplot(pos_df,
                 aes(x = fct_reorder(residue_pos, n, .desc = TRUE), y = n)) +
      geom_col(fill = viridis(1, begin = 0.4)) +
      geom_text(aes(label = sprintf("%d\n(%.0f%%)", n, pct)),
                vjust = -0.3, size = 3) +
      facet_wrap(~subtype, scales = "free_x") +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(
        title    = "Ridge NS substitutions by position",
        subtitle = "Number of independent terminal-branch events per position",
        x = "Mature residue position", y = "Count"
      ) +
      theme_bw(base_size = 12) +
      theme(panel.grid.major.x = element_blank())

    ggsave(file.path(out_dir, "mut_positions.png"),
           p4, width = 9, height = 5, dpi = 180)
    message("  Saved: mut_positions.png")
  }

  # ── 5. AA change convergence bar chart ────────────────────────────────────
  change_df <- mut_rec |>
    count(subtype, change) |>
    group_by(subtype) |>
    mutate(change = fct_reorder(change, n))

  p5 <- ggplot(change_df, aes(x = n, y = change)) +
    geom_col(fill = viridis(1, begin = 0.4)) +
    geom_text(aes(label = n), hjust = -0.2, size = 2.5) +
    facet_wrap(~subtype, scales = "free") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      title    = paste0("Independent occurrences of each ", region, " AA change"),
      subtitle = paste0("Mutations occurring > ", min_count, " times"),
      x = "Count", y = "AA change"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.major.y = element_blank())

  ggsave(file.path(out_dir, "mut_changes.png"),
         p5, width = 11, height = 9, dpi = 180)
  message("  Saved: mut_changes.png")

  # ── 6. Temporal distribution ──────────────────────────────────────────────
  temp_df <- mut_rec |>
    filter(!is.na(year)) |>
    count(subtype, change, year) |>
    group_by(subtype) |>
    mutate(change = fct_reorder(change, n, .fun = sum, .desc = TRUE))

  n_years  <- length(unique(temp_df$year))
  year_pal <- setNames(viridis(n_years, end = 0.9),
                       sort(unique(temp_df$year)))

  p6 <- ggplot(temp_df, aes(x = change, y = n, fill = factor(year))) +
    geom_col(position = "stack") +
    facet_wrap(~subtype, scales = "free_x") +
    scale_fill_manual(values = year_pal, name = "Year") +
    labs(
      title    = paste0("Temporal distribution of recurrent ", region, " mutations"),
      subtitle = paste0("Mutations occurring > ", min_count, " times independently"),
      x = NULL, y = "Count"
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position    = "bottom",
          axis.text.x        = element_text(angle = 45, hjust = 1),
          panel.grid.major.x = element_blank(),
          panel.grid.minor   = element_blank())

  ggsave(file.path(out_dir, "mut_temporal.png"),
         p6, width = 12, height = 7, dpi = 180)
  message("  Saved: mut_temporal.png")

  # ── 7. Chord diagram: AA → AA ────────────────────────────────────────────
  # Link colour: position (ridge) or antigenic region (epitope)
  chord_df <- if (region == "ridge") {
    mut_rec |> count(subtype, residue_pos, parent_aa, tip_aa)
  } else {
    mut_rec |> count(subtype, epi_region, parent_aa, tip_aa)
  }

  link_key  <- if (region == "ridge") "residue_pos" else "epi_region"
  all_keys  <- sort(unique(as.character(chord_df[[link_key]][
    !is.na(chord_df[[link_key]])])))
  link_pal  <- setNames(viridis(length(all_keys), option = "turbo", end = 0.9),
                        all_keys)
  link_labels <- if (region == "ridge") all_keys else
                 sub("epitope_", "Epitope ", all_keys)

  all_aa   <- sort(unique(c(chord_df$parent_aa, chord_df$tip_aa)))
  aa_pal   <- setNames(viridis(length(all_aa), end = 0.9), all_aa)

  png(file.path(out_dir, "mut_chord.png"),
      width = 3000, height = 1500, res = 180)
  par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))

  for (st in c("H1N1", "H3N2")) {
    df <- chord_df |> filter(subtype == st, !is.na(.data[[link_key]]))
    if (nrow(df) == 0) next

    st_aa    <- sort(unique(c(df$parent_aa, df$tip_aa)))
    grid_col <- aa_pal[st_aa]
    link_col <- add_transparency(
      link_pal[as.character(df[[link_key]])], 0.2
    )

    circos.clear()
    circos.par(gap.after = 4, start.degree = 90)
    chordDiagram(
      df[, c("parent_aa", "tip_aa", "n")],
      grid.col = grid_col, col = link_col, transparency = 0,
      directional = 1, direction.type = "arrows", link.arr.type = "big.arrow",
      annotationTrack = "grid", preAllocateTracks = 1
    )
    circos.trackPlotRegion(track.index = 1, bg.border = NA,
      panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.3,
                    get.cell.meta.data("sector.index"),
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = 0.9)
      })
    title(st, cex.main = 1.2)
    legend_title <- if (region == "ridge") "Position" else "Antigenic region"
    legend("bottomleft", legend = link_labels, fill = link_pal,
           title = legend_title, cex = 0.7, bty = "n", ncol = 2)
    circos.clear()
  }
  dev.off()
  message("  Saved: mut_chord.png")

  # ── 8. Chord diagram: nucleotide changes ─────────────────────────────────
  nt_df <- mut_rec |>
    filter(!is.na(parent_codon), !is.na(tip_codon),
           !is.na(.data[[link_key]])) |>
    rowwise() |>
    mutate(nt_changes = list({
      p <- strsplit(parent_codon, "")[[1]]
      t <- strsplit(tip_codon,    "")[[1]]
      idx <- which(p != t)
      if (length(idx) == 0) NULL
      else data.frame(from = p[idx], to = t[idx], stringsAsFactors = FALSE)
    })) |>
    ungroup() |>
    filter(!sapply(nt_changes, is.null)) |>
    tidyr::unnest(nt_changes) |>
    count(subtype, !!sym(link_key), from, to)

  all_nt <- c("A", "C", "G", "T")
  nt_pal <- setNames(viridis(4, end = 0.9), all_nt)

  png(file.path(out_dir, "mut_chord_nt.png"),
      width = 3000, height = 1500, res = 180)
  par(mfrow = c(1, 2), mar = c(1, 1, 2, 1))

  for (st in c("H1N1", "H3N2")) {
    df <- nt_df |> filter(subtype == st)
    if (nrow(df) == 0) next

    link_col <- add_transparency(
      link_pal[as.character(df[[link_key]])], 0.2
    )
    circos.clear()
    circos.par(gap.after = 8, start.degree = 90)
    chordDiagram(
      df[, c("from", "to", "n")],
      grid.col = nt_pal[all_nt], col = link_col, transparency = 0,
      directional = 1, direction.type = "arrows", link.arr.type = "big.arrow",
      order = all_nt, annotationTrack = "grid", preAllocateTracks = 1
    )
    circos.trackPlotRegion(track.index = 1, bg.border = NA,
      panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.3,
                    get.cell.meta.data("sector.index"),
                    facing = "clockwise", niceFacing = TRUE,
                    adj = c(0, 0.5), cex = 1.1, font = 2)
      })
    title(st, cex.main = 1.2)
    legend("bottomleft", legend = link_labels, fill = link_pal,
           title = legend_title, cex = 0.7, bty = "n", ncol = 2)
    circos.clear()
  }
  dev.off()
  message("  Saved: mut_chord_nt.png")
}

message("\nDone.")
