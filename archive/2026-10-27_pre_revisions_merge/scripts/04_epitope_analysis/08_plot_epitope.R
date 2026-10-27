# scripts/04_epitope_analysis/08_plot_epitope.R
#
# PURPOSE: Visualise phylolm results for epitope vs non-epitope dN/dS rates.
#
# INPUT:
#   data/rates/all_phylolm_epitope_results.csv
#   data/rates/{gene}_rates_epitope.csv
#
# OUTPUT:  figures/epitope/{phylolm_forest_epitope,omega_density,omega_vs_temp}.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
})

source(here("config.R"))

out_dir <- here("figures", "epitope")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

res_path <- file.path(PATHS$rates, "all_phylolm_epitope_results.csv")
if (!file.exists(res_path)) stop("Results not found: ", res_path)
res <- read_csv(res_path, show_col_types = FALSE)

predictor_labels <- c(
  scaled_temp             = "Temperature",
  scaled_humidity         = "Relative humidity",
  scaled_delta_temp       = "ΔTemperature",
  scaled_delta_humidity   = "ΔRelative humidity",
  scaled_sampling_density = "Sampling density",
  scaled_pop_density      = "Population density"
)

# Clean response labels: strip log_ prefix, replace _ with space,
# annotate site class
response_labels <- c(
  log_epitope_dN_rate      = "Epitope dN rate",
  log_epitope_dS_rate      = "Epitope dS rate",
  log_epitope_omega        = "Epitope ω",
  log_nonepitope_dN_rate   = "Non-epitope dN rate",
  log_nonepitope_dS_rate   = "Non-epitope dS rate",
  log_nonepitope_omega     = "Non-epitope ω",
  log_ridge_dN_rate        = "Ridge dN rate",
  log_ridge_dS_rate        = "Ridge dS rate",
  log_ridge_omega          = "Ridge ω",
  log_nonridge_dN_rate     = "Non-ridge dN rate",
  log_nonridge_dS_rate     = "Non-ridge dS rate",
  log_nonridge_omega       = "Non-ridge ω"
)

site_class_labels <- c(
  epitope    = "Antigenic sites",
  nonepitope = "Non-epitope sites",
  ridge      = "Antigenic ridge",
  nonridge   = "Non-ridge sites"
)

pred_colours <- setNames(
  viridisLite::viridis(length(predictor_labels), end = 0.92),
  predictor_labels
)

# ── 1. Forest plots: epitope vs non-epitope coefficients ─────────────────────
rate_type_labels <- c(
  dN_rate = "dN rate",
  dS_rate = "dS rate",
  omega   = "Branch ω (dN/dS)"
)

plot_df <- res |>
  filter(term %in% names(predictor_labels),
         response %in% names(response_labels)) |>
  mutate(
    ci_lo      = estimate - 1.96 * std_error,
    ci_hi      = estimate + 1.96 * std_error,
    predictor  = factor(predictor_labels[term], levels = predictor_labels),
    site_class = factor(
      case_when(
        grepl("^log_ridge_",      response) ~ "Antigenic ridge",
        grepl("^log_nonridge_",   response) ~ "Non-ridge sites",
        grepl("^log_epitope_",    response) ~ "Antigenic sites",
        grepl("^log_nonepitope_", response) ~ "Non-epitope sites"
      ),
      levels = c("Antigenic ridge", "Antigenic sites",
                 "Non-epitope sites", "Non-ridge sites")
    ),
    rate_type  = factor(
      case_when(
        grepl("_dN_rate$", response) ~ "dN_rate",
        grepl("_dS_rate$", response) ~ "dS_rate",
        grepl("_omega$",   response) ~ "omega"
      ),
      levels = names(rate_type_labels),
      labels = rate_type_labels
    ),
    subtype    = ifelse(gene %in% c("H1", "N1"), "H1N1", "H3N2")
  )

make_forest <- function(df, sig_col, subtitle, filename) {
  df <- df |>
    mutate(point_fill = ifelse(.data[[sig_col]],
                               pred_colours[as.character(predictor)], "white"))
  ggplot(df, aes(x = estimate, y = site_class,
                 colour = predictor, group = predictor)) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                   linewidth = 0.7,
                   position = position_dodge(width = 0.7)) +
    geom_point(aes(fill = point_fill),
               shape = 21, size = 3, stroke = 0.8,
               position = position_dodge(width = 0.7)) +
    scale_colour_manual(values = pred_colours, name = NULL) +
    scale_fill_identity() +
    scale_y_discrete(limits = rev) +
    facet_grid(rate_type ~ subtype) +
    labs(
      title    = "Climate predictors of epitope vs non-epitope dN/dS rates",
      subtitle = subtitle,
      x = "Coefficient estimate", y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(strip.text         = element_text(face = "bold", size = 10),
          legend.position    = "bottom",
          panel.grid.minor   = element_blank(),
          panel.grid.major.x = element_blank())
}

plot_df_raw <- plot_df |>
  mutate(sig = !is.na(p.value) & p.value < 0.05)

plot_df_fdr <- plot_df |>
  mutate(sig = !is.na(p.fdr) & p.fdr < 0.05)

p1 <- make_forest(plot_df_raw, "sig",
                  "Coefficients ± 95% CI  |  Filled = p < 0.05",
                  "phylolm_forest_epitope.png")
ggsave(file.path(out_dir, "phylolm_forest_epitope.png"),
       p1, width = 10, height = 10, dpi = 180)
message("Saved: phylolm_forest_epitope.png")

p1_fdr <- make_forest(plot_df_fdr, "sig",
                      "Coefficients ± 95% CI  |  Filled = BH-FDR p < 0.05",
                      "phylolm_forest_epitope_fdr.png")
ggsave(file.path(out_dir, "phylolm_forest_epitope_fdr.png"),
       p1_fdr, width = 10, height = 10, dpi = 180)
message("Saved: phylolm_forest_epitope_fdr.png")

# ── 2. Rate density plots ──────────────────────────────────────────────────────
dat <- lapply(c("H1", "H3"), function(g) {
  path <- file.path(PATHS$rates, paste0(g, "_rates_epitope.csv"))
  if (!file.exists(path)) return(NULL)
  read_csv(path, show_col_types = FALSE) |>
    mutate(gene    = g,
           subtype = ifelse(g == "H1", "H1N1", "H3N2"))
}) |> bind_rows()

site_class_levels <- c("ridge", "epitope", "nonepitope", "nonridge")
site_class_labels <- c("Antigenic ridge", "Antigenic sites",
                       "Non-epitope sites", "Non-ridge sites")

site_colours <- c(
  "Antigenic ridge"   = "#FDE725",
  "Antigenic sites"   = "#35B779",
  "Non-epitope sites" = "#31688E",
  "Non-ridge sites"   = "#440154"
)

make_rate_density <- function(dat, rate_suffix, xlab, title, vline = NULL) {
  long <- dat |>
    select(gene, subtype,
           ridge      = paste0("log_ridge_",      rate_suffix),
           epitope    = paste0("log_epitope_",    rate_suffix),
           nonepitope = paste0("log_nonepitope_", rate_suffix),
           nonridge   = paste0("log_nonridge_",   rate_suffix)) |>
    pivot_longer(c(ridge, epitope, nonepitope, nonridge),
                 names_to = "site_class", values_to = "value") |>
    filter(is.finite(value)) |>
    mutate(site_class = factor(site_class,
                               levels = site_class_levels,
                               labels = site_class_labels))

  p <- ggplot(long, aes(x = value, fill = site_class, colour = site_class)) +
    geom_density(alpha = 0.4, linewidth = 0.7) +
    facet_wrap(~subtype) +
    scale_fill_manual(values   = site_colours, name = NULL) +
    scale_colour_manual(values = site_colours, name = NULL) +
    labs(title = title, x = xlab, y = "Density") +
    theme_bw(base_size = 12) +
    theme(legend.position = "bottom")

  if (!is.null(vline))
    p <- p + geom_vline(xintercept = vline, linetype = "dashed",
                        colour = "grey40")
  p
}

p2 <- make_rate_density(dat, "omega",
  xlab  = "log(branch ω)",
  title = "Distribution of branch ω (dN/dS) by site class",
  vline = 0)
ggsave(file.path(out_dir, "omega_density.png"),
       p2, width = 9, height = 5, dpi = 180)
message("Saved: omega_density.png")

p2b <- make_rate_density(dat, "dN_rate",
  xlab  = "log(dN rate)",
  title = "Distribution of branch dN rate by site class")
ggsave(file.path(out_dir, "dN_density.png"),
       p2b, width = 9, height = 5, dpi = 180)
message("Saved: dN_density.png")

p2c <- make_rate_density(dat, "dS_rate",
  xlab  = "log(dS rate)",
  title = "Distribution of branch dS rate by site class")
ggsave(file.path(out_dir, "dS_density.png"),
       p2c, width = 9, height = 5, dpi = 180)
message("Saved: dS_density.png")

# ── 3. Epitope omega vs temperature ───────────────────────────────────────────
gene_colours <- c(H1 = "#E41A1C", H3 = "#377EB8")

omega_vs_temp_long <- dat |>
  select(gene, subtype, era5_temp_mean,
         epitope = log_epitope_omega,
         ridge   = log_ridge_omega) |>
  pivot_longer(c(epitope, ridge),
               names_to = "site_class", values_to = "log_omega") |>
  filter(is.finite(log_omega), is.finite(era5_temp_mean)) |>
  mutate(site_class = factor(site_class,
                             levels = c("ridge", "epitope"),
                             labels = c("Antigenic ridge",
                                        "Antigenic sites")))

p3 <- ggplot(omega_vs_temp_long,
             aes(x = era5_temp_mean, y = log_omega,
                 colour = gene)) +
  geom_point(alpha = 0.1, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE,
              linewidth = 1, alpha = 0.15) +
  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey40") +
  facet_grid(site_class ~ subtype) +
  scale_colour_manual(values = gene_colours, name = NULL) +
  labs(
    title = "Antigenic site branch ω vs temperature",
    x     = "Monthly mean temperature (°C)",
    y     = "log(branch ω)"
  ) +
  theme_bw(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(file.path(out_dir, "omega_vs_temp.png"),
       p3, width = 9, height = 5, dpi = 180)
message("Saved: omega_vs_temp.png")

message("\nDone.")
