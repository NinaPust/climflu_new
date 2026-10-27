# scripts/05_final_figs/fig_jk_a.R
#
# Fig JK-A: Jackknife coefficient distributions
#   Density plot of temperature coefficient estimates across 1,000 jackknife
#   subsamples, H1 (H1N1) and H3 (H3N2) side by side.
#   Vertical line = full-dataset estimate. Rug marks below density.
#   Shows the signal is not driven by a few lucky subsamples.
#
# INPUT:  data/rates/jackknife/{H1,H3}/jackknife_results.csv
#         data/rates/phylolm/all_phylolm_results.csv  (full-dataset estimate)
# OUTPUT: paper/figures/fig_jk_a.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
})

source(here("config.R"))

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load jackknife results ────────────────────────────────────────────────────
jk <- bind_rows(
  read_csv(file.path(PATHS$rates, "jackknife", "H1",
                     "jackknife_results.csv"), show_col_types = FALSE),
  read_csv(file.path(PATHS$rates, "jackknife", "H3",
                     "jackknife_results.csv"), show_col_types = FALSE)
) |>
  mutate(
    subtype = ifelse(gene == "H1", "H1N1 (H1)", "H3N2 (H3)"),
    sig     = case_when(
      p.value < 0.001 ~ "p < 0.001",
      p.value < 0.01  ~ "p < 0.01",
      p.value < 0.05  ~ "p < 0.05",
      TRUE            ~ "ns"
    )
  )

# ── Full-dataset estimate ─────────────────────────────────────────────────────
full_est <- read_csv(
  file.path(PATHS$rates_phylolm, "all_phylolm_results.csv"),
  show_col_types = FALSE
) |>
  filter(response == "log_rate", term == "scaled_temp",
         gene %in% c("H1", "H3")) |>
  mutate(subtype = ifelse(gene == "H1", "H1N1 (H1)", "H3N2 (H3)"))

n_samples <- n_distinct(jk$sample_id)
n_per     <- unique(jk$n)[1]

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(jk, aes(x = estimate)) +

  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +

  geom_vline(data = full_est,
             aes(xintercept = estimate),
             colour = viridis(3)[2], linewidth = 0.9, linetype = "solid") +

  geom_density(fill = viridis(3)[1], colour = viridis(3)[1],
               alpha = 0.35, linewidth = 0.7) +

  geom_rug(alpha = 0.15, linewidth = 0.3) +

  facet_wrap(~ subtype, nrow = 1, scales = "free_y") +

  labs(
    x       = "Temperature coefficient (jackknife subsamples)",
    y       = "Density",
    caption = paste0(
      n_samples, " subsamples of n = ", n_per,
      " sequences each | Vertical line = full-dataset estimate"
    )
  ) +

  theme_bw(base_size = 13) +
  theme(
    strip.text         = element_text(face = "bold", size = 14),
    strip.background   = element_blank(),
    panel.grid.minor   = element_blank(),
    axis.title         = element_text(size = 14),
    axis.text          = element_text(size = 12),
    aspect.ratio       = 1,
    plot.caption       = element_text(size = 9, colour = "grey40")
  )

ggsave(file.path(out_dir, "fig_jk_a.png"),
       p, width = 9, height = 5.5, dpi = 200)
message("Saved: paper/figures/fig_jk_a.png")
