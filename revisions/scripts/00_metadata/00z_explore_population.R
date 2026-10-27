# scripts/00_metadata/00z_explore_population.R
#
# PURPOSE: Exploratory diagnostics on GHS-POP population density data.
#          Gitignored — not part of the reproducible pipeline.
#
# INPUT:  data/population/{gene}_popdens.csv  (run 03_population_density.py first)
# OUTPUT: figures/population_diagnostics/     (PNG / PDF files)
#
# USAGE:  Rscript scripts/00_metadata/00z_explore_population.R

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

out_dir <- here("revisions", "figures", "population_diagnostics")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

URBAN_THRESHOLD <- 500   # people/km²; see CLAUDE.md
ALL_GENES <- c(GENES, "PA_H1N1", "PA_H3N2")

# ── Load all population density data ─────────────────────────────────────────
pop_list <- lapply(ALL_GENES, function(g) {
  path <- file.path(PATHS$population, paste0(g, "_popdens.csv"))
  if (!file.exists(path)) {
    warning("Missing: ", path)
    return(NULL)
  }
  read_csv(path, show_col_types = FALSE) |>
    select(accession, latitude, longitude, year, month, country, region,
           pop_density, pop_epoch) |>
    mutate(gene = g)
})
pop <- bind_rows(pop_list)

if (nrow(pop) == 0) stop("No data loaded — run 03_population_density.py first.")

pop <- pop |>
  mutate(
    urban_rural  = case_when(
      is.na(pop_density)            ~ NA_character_,
      pop_density >= URBAN_THRESHOLD ~ "Urban",
      TRUE                           ~ "Rural"
    ),
    log10_density = log10(pop_density + 1)
  )

message("Loaded ", nrow(pop), " records across ", n_distinct(pop$gene), " genes")
message("  NA pop_density: ", sum(is.na(pop$pop_density)),
        " (", round(100 * mean(is.na(pop$pop_density)), 1), "%)")

# ── 1. Missing data summary ───────────────────────────────────────────────────
na_summary <- pop |>
  group_by(gene) |>
  summarise(
    n_total   = n(),
    n_na      = sum(is.na(pop_density)),
    pct_na    = round(100 * mean(is.na(pop_density)), 1),
    pct_urban = round(100 * mean(urban_rural == "Urban", na.rm = TRUE), 1),
    pct_rural = round(100 * mean(urban_rural == "Rural", na.rm = TRUE), 1),
    .groups = "drop"
  )
message("\n── Missing / urban-rural summary ─────────────────────────────")
print(na_summary, n = Inf)
write_csv(na_summary, file.path(out_dir, "na_summary.csv"))

# ── 2. Density distribution (log scale) ──────────────────────────────────────
p_hist <- pop |>
  filter(!is.na(pop_density), pop_density > 0) |>
  ggplot(aes(x = log10_density, fill = gene)) +
  geom_histogram(bins = 60, alpha = 0.7, position = "identity") +
  geom_vline(xintercept = log10(URBAN_THRESHOLD + 1),
             linetype = "dashed", colour = "black", linewidth = 0.6) +
  annotate("text", x = log10(URBAN_THRESHOLD + 1) + 0.05,
           y = Inf, label = "Urban threshold", hjust = 0, vjust = 1.5,
           size = 3) +
  scale_fill_brewer(palette = "Set2") +
  facet_wrap(~gene, scales = "free_y", ncol = 3) +
  labs(title = "Population density distribution by gene segment",
       subtitle = "Sequences with pop_density > 0 (log₁₀ scale)",
       x = "log₁₀(people/km²)", y = "Count", fill = "Segment") +
  theme_minimal(base_size = 11) +
  theme(legend.position = "none")

ggsave(file.path(out_dir, "01_density_distribution.png"), p_hist,
       width = 11, height = 7, dpi = 150)
message("Saved: 01_density_distribution.png")

# ── 3. Urban vs rural breakdown per gene ─────────────────────────────────────
p_urban <- pop |>
  filter(!is.na(urban_rural)) |>
  count(gene, urban_rural) |>
  group_by(gene) |>
  mutate(pct = n / sum(n)) |>
  ungroup() |>
  ggplot(aes(x = gene, y = pct, fill = urban_rural)) +
  geom_col() +
  geom_text(aes(label = paste0(round(pct * 100), "%")),
            position = position_stack(vjust = 0.5), size = 3.5) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_manual(values = c(Urban = "#e07b39", Rural = "#4a9e6b")) +
  labs(title = paste0("Urban/rural split (threshold: >",
                      URBAN_THRESHOLD, " people/km²)"),
       x = "Gene segment", y = "Proportion", fill = NULL) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "top")

ggsave(file.path(out_dir, "02_urban_rural_split.png"), p_urban,
       width = 8, height = 5, dpi = 150)
message("Saved: 02_urban_rural_split.png")

# ── 4. Density by WHO region ──────────────────────────────────────────────────
p_region <- pop |>
  filter(!is.na(pop_density), !is.na(region), region != "") |>
  mutate(region = fct_reorder(region, pop_density, median, na.rm = TRUE)) |>
  ggplot(aes(x = log10_density, y = region, fill = region)) +
  geom_boxplot(outlier.size = 0.4, outlier.alpha = 0.3, show.legend = FALSE) +
  geom_vline(xintercept = log10(URBAN_THRESHOLD + 1),
             linetype = "dashed", colour = "grey30", linewidth = 0.5) +
  scale_fill_brewer(palette = "Set3") +
  labs(title = "Population density of sample locations by WHO region",
       x = "log₁₀(people/km²)", y = NULL) +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "03_density_by_region.png"), p_region,
       width = 9, height = 5, dpi = 150)
message("Saved: 03_density_by_region.png")

# ── 5. Density over time (median per year) ────────────────────────────────────
p_time <- pop |>
  filter(!is.na(pop_density), !is.na(year)) |>
  group_by(gene, year) |>
  summarise(median_dens = median(pop_density, na.rm = TRUE), .groups = "drop") |>
  ggplot(aes(x = year, y = log10(median_dens + 1), colour = gene, group = gene)) +
  geom_line() +
  geom_point(size = 1.5) +
  scale_colour_brewer(palette = "Set2") +
  labs(title = "Median sample population density over time",
       x = "Year", y = "log₁₀(median people/km²)", colour = "Segment") +
  theme_minimal(base_size = 12)

ggsave(file.path(out_dir, "04_density_over_time.png"), p_time,
       width = 9, height = 5, dpi = 150)
message("Saved: 04_density_over_time.png")

# ── 6. Top countries by median population density ────────────────────────────
top_countries <- pop |>
  filter(!is.na(pop_density), !is.na(country)) |>
  group_by(country) |>
  summarise(n = n(), median_dens = median(pop_density, na.rm = TRUE),
            .groups = "drop") |>
  filter(n >= 30) |>
  slice_max(median_dens, n = 20)

p_country_dens <- top_countries |>
  mutate(country = fct_reorder(country, median_dens)) |>
  ggplot(aes(x = median_dens, y = country)) +
  geom_col(fill = "#5b8fcc") +
  geom_text(aes(label = paste0("n=", n)), hjust = -0.1, size = 3) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
  labs(title = "Top 20 countries by median sample population density",
       subtitle = "Countries with ≥30 sequences; n = number of sequences",
       x = "Median people/km²", y = NULL) +
  theme_minimal(base_size = 11)

ggsave(file.path(out_dir, "05_top_countries_density.png"), p_country_dens,
       width = 9, height = 7, dpi = 150)
message("Saved: 05_top_countries_density.png")

# ── 7. World map coloured by population density ───────────────────────────────
message("\nGenerating sample maps...")
world <- ne_countries(scale = "medium", returnclass = "sf")

make_density_map <- function(df, gene_label) {
  df <- df |> filter(!is.na(latitude), !is.na(longitude), !is.na(pop_density))
  ggplot() +
    geom_sf(data = world, fill = "grey92", colour = "grey70", linewidth = 0.1) +
    geom_point(data = df |> arrange(pop_density),
               aes(x = longitude, y = latitude, colour = log10(pop_density + 1)),
               size = 0.5, alpha = 0.6) +
    scale_colour_viridis_c(
      option  = "inferno",
      name    = "log₁₀\n(people/km²)",
      limits  = c(0, 4.5),
      oob     = scales::squish
    ) +
    labs(title = paste(gene_label, "— sample locations coloured by pop density")) +
    theme_minimal(base_size = 10) +
    theme(legend.position  = "right",
          panel.grid       = element_blank(),
          axis.text        = element_blank(),
          axis.title       = element_blank())
}

pdf(file.path(out_dir, "06_sample_maps_density.pdf"), width = 11, height = 6)
for (g in ALL_GENES) {
  df_g <- pop |> filter(gene == g)
  print(make_density_map(df_g, g))
}
dev.off()
message("Saved: 06_sample_maps_density.pdf")

# ── 8. Combined summary panel (H1 example) ───────────────────────────────────
h1 <- pop |> filter(gene == "H1", !is.na(pop_density))

pa_hist <- h1 |>
  filter(pop_density > 0) |>
  ggplot(aes(x = log10_density, fill = urban_rural)) +
  geom_histogram(bins = 50, alpha = 0.85) +
  scale_fill_manual(values = c(Urban = "#e07b39", Rural = "#4a9e6b"),
                    na.value = "grey70") +
  geom_vline(xintercept = log10(URBAN_THRESHOLD + 1),
             linetype = "dashed", linewidth = 0.6) +
  labs(x = "log₁₀(people/km²)", y = "Count", fill = NULL,
       title = "Density distribution") +
  theme_minimal(base_size = 11)

pa_region <- h1 |>
  filter(!is.na(region)) |>
  mutate(region = fct_reorder(region, pop_density, median)) |>
  ggplot(aes(x = log10_density, y = region, fill = urban_rural)) +
  geom_boxplot(outlier.size = 0.5, outlier.alpha = 0.3, show.legend = FALSE) +
  scale_fill_manual(values = c(Urban = "#e07b39", Rural = "#4a9e6b"),
                    na.value = "grey70") +
  labs(x = "log₁₀(people/km²)", y = NULL, title = "By WHO region") +
  theme_minimal(base_size = 11)

p_summary <- pa_hist + pa_region +
  plot_annotation(title = "H1 — population density of sample locations",
                  theme = theme(plot.title = element_text(size = 14, face = "bold")))

ggsave(file.path(out_dir, "07_H1_summary_panel.png"), p_summary,
       width = 12, height = 5, dpi = 150)
message("Saved: 07_H1_summary_panel.png")

message("\nAll diagnostics written to: ", out_dir)
