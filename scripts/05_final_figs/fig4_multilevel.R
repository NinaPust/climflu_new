# scripts/05_final_figs/fig4_multilevel.R
#
# Fig 4: Temperature effect is consistent across genomic regions
#   Forest plot of scaled_temp coefficient (± 95% CI) for five response
#   variables, ordered from whole-timetree rate down to the antigenic ridge:
#     1. Timetree log(branch rate)         [all_phylolm_results.csv]
#     2. Full-gene log(NS rate)             [all_phylolm_ns_results.csv]
#     3. Full-gene log(S rate)              [all_phylolm_ns_results.csv]
#     4. Epitope log(NS rate)              [all_phylolm_ns_results.csv]
#     5. Ridge log(NS rate)                [all_phylolm_ns_results.csv]
#   H1N1 and H3N2 side by side.
#   Filled points = BH-FDR p < 0.05.
#
# INPUT:  data/rates/phylolm/all_phylolm_results.csv
#         data/rates/phylolm/all_phylolm_ns_results.csv
# OUTPUT: paper/figures/fig4_multilevel.png

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

# ── Load ──────────────────────────────────────────────────────────────────────
main_res <- read_csv(file.path(PATHS$rates_phylolm, "all_phylolm_results.csv"),
                     show_col_types = FALSE) |>
  filter(response == "log_rate", term == "scaled_temp",
         gene %in% c("H1", "H3")) |>
  mutate(region = "Terminal\nbranch rate")

ns_res <- read_csv(file.path(PATHS$rates_phylolm, "all_phylolm_ns_results.csv"),
                   show_col_types = FALSE) |>
  filter(term == "scaled_temp", gene %in% c("H1", "H3")) |>
  mutate(region = case_when(
    response == "log_full_ns_rate"    ~ "Whole-gene\nNS rate",
    response == "log_full_s_rate"     ~ "Whole-gene\nS rate",
    response == "log_epitope_ns_rate" ~ "Epitope\nNS rate",
    response == "log_ridge_ns_rate"   ~ "Ridge\nNS rate",
    TRUE                              ~ NA_character_
  )) |>
  filter(!is.na(region))

# Combine
combined <- bind_rows(main_res, ns_res) |>
  mutate(
    region  = factor(region, levels = c(
      "Terminal\nbranch rate",
      "Whole-gene\nNS rate",
      "Whole-gene\nS rate",
      "Epitope\nNS rate",
      "Ridge\nNS rate"
    )),
    subtype   = ifelse(gene == "H1", "H1N1", "H3N2"),
    ci_lo     = estimate - 1.96 * std_error,
    ci_hi     = estimate + 1.96 * std_error,
    sig_fdr   = !is.na(p.fdr) & p.fdr < 0.05,
    sig_nom   = p.value < 0.05,
    # colour by region type: rate vs NS vs S
    rate_type = case_when(
      region == "Terminal\nbranch rate" ~ "Branch rate",
      grepl("NS", region)               ~ "NS rate",
      grepl("S rate", region)           ~ "S rate"
    ),
    rate_type = factor(rate_type, levels = c("Branch rate", "NS rate", "S rate"))
  )

TYPE_COLOURS <- setNames(
  viridis(3, end = 0.85),
  c("Branch rate", "NS rate", "S rate")
)

p <- ggplot(combined,
            aes(x = estimate, y = region, colour = rate_type)) +


  geom_vline(xintercept = 0, linetype = "dashed",
             colour = "grey50", linewidth = 0.4) +

  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi),
                 linewidth = 1) +

  geom_point(aes(fill = ifelse(sig_fdr, TYPE_COLOURS[as.character(rate_type)], "white")),
             shape = 21, size = 4, stroke = 1) +

  scale_colour_manual(values = TYPE_COLOURS, name = "Rate type") +
  scale_fill_identity() +
  scale_y_discrete(limits = rev) +

  facet_wrap(~ subtype, nrow = 1) +

  labs(x = "Coefficient", y = NULL) +

  theme_bw(base_size = 12) +
  theme(
    strip.text         = element_text(face = "bold", size = 14),
    strip.background   = element_blank(),
    legend.position    = "bottom",
    legend.text        = element_text(size = 12),
    panel.grid.minor   = element_blank(),
    panel.grid.major.x = element_blank(),
    axis.title.x       = element_text(size = 14),
    axis.text          = element_text(size = 12),
    aspect.ratio       = 1
  )

ggsave(file.path(out_dir, "fig4_multilevel.png"),
       p, width = 9, height = 6, dpi = 200)
message("Saved: paper/figures/fig4_multilevel.png")
