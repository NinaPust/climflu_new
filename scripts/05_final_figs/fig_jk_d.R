# scripts/05_final_figs/fig_jk_d.R
#
# Fig JK-D: Significance proportion summary
#   Stacked bar showing what fraction of 1,000 jackknife subsamples falls
#   into each significance tier (ns / p<0.05 / p<0.01 / p<0.001),
#   for H1 (H1N1) and H3 (H3N2).
#   Communicates overall robustness at a glance.
#
# INPUT:  data/rates/jackknife/{H1,H3}/jackknife_results.csv
# OUTPUT: paper/figures/fig_jk_d.png

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(tidyr)
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

# ── Load and summarise ────────────────────────────────────────────────────────
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
  )

prop_df <- jk |>
  count(subtype, sig) |>
  group_by(subtype) |>
  mutate(prop = n / sum(n)) |>
  ungroup()

n_total <- n_distinct(jk$sample_id)

# ── Plot ──────────────────────────────────────────────────────────────────────
p <- ggplot(prop_df, aes(x = subtype, y = prop, fill = sig)) +

  geom_col(width = 0.55, colour = "white", linewidth = 0.3) +

  geom_text(aes(label = ifelse(prop > 0.04,
                               paste0(round(prop * 100, 1), "%"), "")),
            position = position_stack(vjust = 0.5),
            size = 4, colour = "white", fontface = "bold") +

  scale_fill_manual(values = SIG_COLOURS, name = NULL,
                    guide = guide_legend(reverse = TRUE)) +
  scale_y_continuous(labels = scales::percent_format(),
                     expand = expansion(mult = c(0, 0.02))) +

  labs(
    x       = NULL,
    y       = paste0("Proportion of subsamples (n = ", n_total, ")"),
    caption = "Each subsample: 1,000 sequences drawn without replacement"
  ) +

  theme_bw(base_size = 13) +
  theme(
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    legend.position    = "right",
    legend.text        = element_text(size = 12),
    axis.title.y       = element_text(size = 13),
    axis.text          = element_text(size = 13),
    aspect.ratio       = 1.4,
    plot.caption       = element_text(size = 9, colour = "grey40")
  )

ggsave(file.path(out_dir, "fig_jk_d.png"),
       p, width = 6, height = 6, dpi = 200)
message("Saved: paper/figures/fig_jk_d.png")
