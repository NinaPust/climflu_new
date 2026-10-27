# scripts/05_final_figs/fig_jk_b.R
#
# Fig JK-B: Jackknife ordered strip plot
#   All 1,000 subsamples sorted by coefficient estimate (ascending), with
#   error bars (±1 SE), coloured by significance tier.
#   H1N1 and H3N2 stacked vertically.
#
# INPUT:  data/rates/jackknife/{H1,H3}/jackknife_results.csv
# OUTPUT: paper/figures/fig_jk_b.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridis)
})

source(here("config.R"))
source(here("R", "jackknife_utils.R"))

# Coefficient to pull from the long-format jackknife results
# (one row per replicate per term). Override per run if needed.
JK_COEF <- "scaled_temp"


out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load and prepare ──────────────────────────────────────────────────────────
SIG_LEVELS  <- c("ns", "p < 0.05", "p < 0.01", "p < 0.001")
SIG_COLOURS <- setNames(viridis(4, end = 0.92), SIG_LEVELS)

jk <- bind_rows(
  read_jackknife(file.path(PATHS$rates, "jackknife", "H1",
                           "jackknife_results.csv"), JK_COEF),
  read_jackknife(file.path(PATHS$rates, "jackknife", "H3",
                           "jackknife_results.csv"), JK_COEF)
) |>
  mutate(
    subtype = ifelse(gene == "H1", "H1N1 (H1)", "H3N2 (H3)"),
    sig     = factor(
      case_when(
        p.value < 0.001 ~ "p < 0.001",
        p.value < 0.01  ~ "p < 0.01",
        p.value < 0.05  ~ "p < 0.05",
        TRUE            ~ "ns"
      ),
      levels = SIG_LEVELS
    )
  ) |>
  group_by(subtype) |>
  arrange(estimate) |>
  mutate(rank = row_number()) |>
  ungroup()

# Proportion significant per subtype (for subtitle)
sig_summary <- jk |>
  group_by(subtype) |>
  summarise(
    pct_sig = round(100 * mean(p.value < 0.05), 1),
    .groups = "drop"
  )

subtitle_text <- paste(
  paste0(sig_summary$subtype, ": ", sig_summary$pct_sig, "% p < 0.05"),
  collapse = "   |   "
)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(jk, aes(x = rank, y = estimate, colour = sig)) +

  geom_hline(yintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +

  geom_errorbar(aes(ymin = estimate - std_error,
                    ymax = estimate + std_error),
                width = 0, alpha = 0.4, linewidth = 0.25) +

  geom_point(size = 0.8, alpha = 0.8) +

  scale_colour_manual(values = SIG_COLOURS, name = NULL,
                      guide = guide_legend(
                        override.aes = list(size = 3, alpha = 1)
                      )) +

  facet_wrap(~ subtype, ncol = 1, scales = "free_x") +

  labs(
    x        = "Sample rank (ordered by estimate)",
    y        = "Temperature coefficient",
    subtitle = subtitle_text
  ) +

  theme_bw(base_size = 13) +
  theme(
    strip.text       = element_text(face = "bold", size = 14),
    strip.background = element_blank(),
    legend.position  = "bottom",
    legend.text      = element_text(size = 12),
    panel.grid.minor = element_blank(),
    axis.title       = element_text(size = 14),
    axis.text        = element_text(size = 11),
    plot.subtitle    = element_text(size = 10, colour = "grey30")
  )

ggsave(file.path(out_dir, "fig_jk_b.png"),
       p, width = 10, height = 7, dpi = 200)
message("Saved: paper/figures/fig_jk_b.png")
