# scripts/00_metadata/00z_explore_metadata.R
#
# PURPOSE: Exploratory diagnostics on parsed GISAID metadata CSVs.
#          Gitignored — not part of the reproducible pipeline.
#
# INPUT:  data/metadata/{gene}_metadata.csv  (run 00_parse_gisaid.R first)
# OUTPUT: figures/metadata_diagnostics/      (PNG files)
#
# USAGE:  Rscript scripts/00_metadata/00z_explore_metadata.R

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(forcats)
  library(patchwork)
  library(sf)
  library(rnaturalearth)
})

source(here("revisions", "config.R"))

out_dir <- here("revisions", "figures", "metadata_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

genes_plot <- c("H1", "H3", "N1", "N2", "PA_H1N1", "PA_H3N2")

# ── Load all metadata ─────────────────────────────────────────────────────────
meta_list <- lapply(genes_plot, function(g) {
  path <- file.path(PATHS$metadata, paste0(g, "_metadata.csv"))
  if (!file.exists(path)) {
    warning("Missing: ", path)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE) |> mutate(gene = g)
})
meta <- bind_rows(meta_list)

if (nrow(meta) == 0) stop("No metadata loaded — run 00_parse_gisaid.R first.")

gene_counts <- paste(genes_plot, sapply(meta_list, nrow),
                     sep = "=", collapse = ", ")
message("Total records loaded: ", nrow(meta), "  (", gene_counts, ")")

# ── 1. Sequences per year ─────────────────────────────────────────────────────
p_year <- meta |>
  filter(!is.na(year)) |>
  ggplot(aes(x = year, fill = gene)) +
  geom_bar(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Sequences per year by gene segment",
       x = "Year", y = "Count", fill = "Segment") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "01_sequences_per_year.png"), p_year,
       width = 10, height = 5, dpi = 150)
message("Saved: 01_sequences_per_year.png")

# ── 2. Top countries ──────────────────────────────────────────────────────────
top_n <- 25

country_counts <- meta |>
  filter(!is.na(country), country != "") |>
  count(gene, country) |>
  group_by(gene) |>
  mutate(rank = rank(-n, ties.method = "first")) |>
  ungroup()

top_countries <- country_counts |>
  group_by(country) |>
  summarise(total = sum(n)) |>
  slice_max(total, n = top_n) |>
  pull(country)

p_country <- country_counts |>
  filter(country %in% top_countries) |>
  mutate(country = fct_reorder(country, n, sum)) |>
  ggplot(aes(x = n, y = country, fill = gene)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = paste("Top", top_n, "countries by sequence count"),
       x = "Count", y = NULL, fill = "Segment") +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir, "02_top_countries.png"), p_country,
       width = 10, height = 8, dpi = 150)
message("Saved: 02_top_countries.png")

# ── 3. WHO region breakdown ───────────────────────────────────────────────────
p_region <- meta |>
  filter(!is.na(region), region != "") |>
  count(gene, region) |>
  mutate(region = fct_reorder(region, n, sum)) |>
  ggplot(aes(x = n, y = region, fill = gene)) +
  geom_col(position = "dodge") +
  scale_fill_brewer(palette = "Set2") +
  labs(title = "Sequences per WHO region",
       x = "Count", y = NULL, fill = "Segment") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "03_region_breakdown.png"), p_region,
       width = 9, height = 5, dpi = 150)
message("Saved: 03_region_breakdown.png")

# ── 4. Seasonal pattern (month) ───────────────────────────────────────────────
p_month <- meta |>
  filter(!is.na(month)) |>
  count(gene, month) |>
  mutate(month = factor(month, levels = 1:12,
                        labels = month.abb)) |>
  ggplot(aes(x = month, y = n, colour = gene, group = gene)) +
  geom_line() + geom_point() +
  scale_colour_brewer(palette = "Set2") +
  labs(title = "Seasonal distribution (all years combined)",
       x = "Month", y = "Count", colour = "Segment") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "04_seasonal_pattern.png"), p_month,
       width = 9, height = 5, dpi = 150)
message("Saved: 04_seasonal_pattern.png")

# ── 5. Missing data summary ───────────────────────────────────────────────────
missing_summary <- meta |>
  group_by(gene) |>
  summarise(
    n_total          = n(),
    pct_no_country   = round(100 * mean(is.na(country) | country == ""), 1),
    pct_no_date      = round(100 * mean(is.na(collection_date)), 1),
    pct_no_month     = round(100 * mean(is.na(month)), 1),
    pct_no_lat       = round(100 * mean(is.na(latitude)), 1),
    pct_no_lineage   = round(100 * mean(is.na(lineage) | lineage == ""), 1),
    .groups = "drop"
  )

message("\n── Missing data summary ──────────────────────────────────────")
print(missing_summary, n = Inf)
write_csv(missing_summary,
          file.path(out_dir, "missing_data_summary.csv"))
message("Saved: missing_data_summary.csv")

# ── 6. Lineage breakdown (H1 only — pdm09 vs seasonal) ───────────────────────
if ("H1" %in% genes_plot) {
  p_lineage <- meta |>
    filter(gene == "H1", !is.na(lineage)) |>
    count(year, lineage) |>
    ggplot(aes(x = year, y = n, fill = lineage)) +
    geom_col() +
    labs(title = "H1 lineage composition over time",
         x = "Year", y = "Count", fill = "Lineage") +
    theme_minimal(base_size = 12)

  ggsave(file.path(out_dir, "05_H1_lineage.png"), p_lineage,
         width = 9, height = 5, dpi = 150)
  message("Saved: 05_H1_lineage.png")
}

# ── 7. Sample maps coloured by collection year ────────────────────────────────
message("\nGenerating sample maps...")

world <- ne_countries(scale = "medium", returnclass = "sf")

make_sample_map <- function(df, gene_label) {
  df <- df |> filter(!is.na(latitude), !is.na(longitude), !is.na(year))
  ggplot() +
    geom_sf(data = world, fill = "grey90", colour = "grey70", linewidth = 0.1) +
    geom_point(data = df,
               aes(x = longitude, y = latitude, colour = year),
               size = 0.5, alpha = 0.5) +
    scale_colour_viridis_c(option = "plasma") +
    labs(title = paste(gene_label, "samples"),
         colour = "Year") +
    theme_minimal(base_size = 10) +
    theme(legend.position = "bottom",
          panel.grid = element_blank())
}

# H3N2 subtype: H3, N2, PA_H3N2
h3_map <- make_sample_map(filter(meta, gene == "H3"), "H3")
n2_map <- make_sample_map(filter(meta, gene == "N2"), "N2")
pa_h3n2_map <- make_sample_map(filter(meta, gene == "PA_H3N2"), "PA (H3N2)")

pdf(file.path(out_dir, "h3n2_map.pdf"), width = 10, height = 6)
print(h3_map)
print(n2_map)
print(pa_h3n2_map)
dev.off()
message("Saved: h3n2_map.pdf")

# H1N1 subtype: H1, N1, PA_H1N1
h1_map <- make_sample_map(filter(meta, gene == "H1"), "H1")
n1_map <- make_sample_map(filter(meta, gene == "N1"), "N1")
pa_h1n1_map <- make_sample_map(filter(meta, gene == "PA_H1N1"), "PA (H1N1)")

pdf(file.path(out_dir, "h1n1_map.pdf"), width = 10, height = 6)
print(h1_map)
print(n1_map)
print(pa_h1n1_map)
dev.off()
message("Saved: h1n1_map.pdf")

message("\nAll diagnostics written to: ", out_dir)
