# scripts/03_rates_analysis/05_brms_multiseed_plots.R
#
# PURPOSE: Diagnostic plots for brms multiseed stability analysis.
#   Page 1 — Forest plot: per-seed posterior means + 95% CrIs, coloured by
#             seed; pooled mean as semi-transparent black diamond.
#   Page 2 — Sample composition: year histogram and lat/lon scatter per seed,
#             to check whether stratified subsampling produced representative
#             draws from the full dataset.
#
# INPUT:  data/rates/brms/{gene}_brms_multiseed.csv
#         data/rates/{gene}_rates.csv   (for sample composition)
# OUTPUT: figures/brms_multiseed_{gene}.pdf
#
# FLAGS:  --gene X   (default: H1)

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
})

source(here("revisions", "config.R"))

args <- commandArgs(trailingOnly = TRUE)
gene <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else "H1"

in_path    <- here("revisions", "data", "rates", "brms",
                   paste0(gene, "_brms_multiseed.csv"))
acc_path   <- here("revisions", "data", "rates", "brms",
                   paste0(gene, "_brms_multiseed_accessions.csv"))
rates_path <- here("revisions", "data", "rates", paste0(gene, "_rates.csv"))
out_path   <- here("revisions", "figures", paste0("brms_multiseed_", gene, ".pdf"))

if (!file.exists(in_path)) stop("File not found: ", in_path)

PREDS <- c("scaled_temp", "scaled_humidity",
           "scaled_delta_temp", "scaled_delta_humidity",
           "scaled_sampling_density", "scaled_pop_density")

PRED_LABELS <- c(
  scaled_temp              = "Temperature",
  scaled_humidity          = "Humidity",
  scaled_delta_temp        = "Delta Temperature",
  scaled_delta_humidity    = "Delta Humidity",
  scaled_sampling_density  = "Sampling density",
  scaled_pop_density       = "Population density"
)

RESP_LABELS <- c(log_rate = "Evolutionary rate (log)", log_omega = "dN/dS (log omega)")

raw <- read_csv(in_path, show_col_types = FALSE) |>
  filter(term %in% PREDS) |>
  mutate(
    term     = factor(term, levels = rev(PREDS), labels = rev(PRED_LABELS[PREDS])),
    response = factor(response, levels = c("log_rate", "log_omega"),
                      labels = RESP_LABELS),
    seed     = factor(seed)
  )

pooled <- raw |>
  group_by(response, term) |>
  summarise(
    mean_pooled = mean(mean),
    sd_pooled   = sqrt(mean(sd^2) + var(mean)),
    .groups = "drop"
  ) |>
  mutate(
    q0.025 = mean_pooled - 1.96 * sd_pooled,
    q0.975 = mean_pooled + 1.96 * sd_pooled
  )

n_seeds <- length(unique(raw$seed))

p <- ggplot() +
  # zero line
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  # per-seed CIs
  geom_linerange(
    data = raw,
    aes(xmin = q0.025, xmax = q0.975, y = term, colour = seed),
    position = position_dodge(width = 0.7),
    linewidth = 0.5, alpha = 0.7
  ) +
  # per-seed means — filled if significant, hollow if not
  geom_point(
    data = raw,
    aes(x = mean, y = term, colour = seed,
        shape = factor(sig_95)),
    position = position_dodge(width = 0.7),
    size = 1.8, alpha = 0.9
  ) +
  # pooled mean diamond
  geom_point(
    data = pooled,
    aes(x = mean_pooled, y = term),
    shape = 18, size = 4, colour = "black", alpha = 0.5
  ) +
  facet_wrap(~response, ncol = 2, scales = "free_x") +
  scale_shape_manual(
    name   = "Significant (95% CrI)",
    values = c("0" = 1, "1" = 19),
    labels = c("0" = "No", "1" = "Yes")
  ) +
  scale_colour_viridis_d(name = "Seed", option = "turbo") +
  labs(
    title    = paste0(gene, " - brms multiseed stability (n ~1000 per seed)"),
    subtitle = paste0(n_seeds, " seeds; black diamond = pooled mean"),
    x        = "Posterior mean (standardised units)",
    y        = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    panel.grid.minor  = element_blank(),
    strip.background  = element_rect(fill = "grey92"),
    legend.position   = "bottom",
    legend.key.size   = unit(0.5, "lines")
  )

# ── Page 2: sample composition per seed ──────────────────────────────────────
p2 <- NULL
if (file.exists(acc_path) && file.exists(rates_path)) {
  acc   <- read_csv(acc_path, show_col_types = FALSE) |>
    mutate(seed = factor(seed))
  rates <- read_csv(rates_path, show_col_types = FALSE)

  samples <- acc |>
    left_join(rates |> select(accession, log_rate, log_omega),
              by = "accession")

  samples <- samples |> mutate(seed = factor(seed))

  p_rate <- ggplot(samples, aes(x = seed, y = log_rate, fill = seed)) +
    geom_violin(alpha = 0.6, colour = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.3,
                 fill = "white", colour = "grey30") +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    labs(title = "Evolutionary rate distribution per seed",
         x = "Seed", y = "log(rate)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  p_omega <- ggplot(
    samples |> filter(is.finite(log_omega)),
    aes(x = seed, y = log_omega, fill = seed)
  ) +
    geom_violin(alpha = 0.6, colour = NA) +
    geom_boxplot(width = 0.15, outlier.size = 0.3,
                 fill = "white", colour = "grey30") +
    scale_fill_viridis_d(option = "turbo", guide = "none") +
    labs(title = "dN/dS distribution per seed",
         x = "Seed", y = "log(omega)") +
    theme_bw(base_size = 10) +
    theme(panel.grid.minor = element_blank())

  p2 <- p_rate + p_omega
} else {
  message("  Accessions file not found, skipping sample composition plot.")
  message("  Re-run 04b_brms_multiseed.R to generate it.")
}

# ── Save ──────────────────────────────────────────────────────────────────────
pdf(out_path, width = 9, height = 5)
print(p)
if (!is.null(p2)) print(p2)
invisible(dev.off())
message("Written: ", out_path)
