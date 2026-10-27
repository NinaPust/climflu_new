# revisions/scripts/05_final_figs/supp9_phyloglm_occurrence.R
#
# Supp 9: Phylogenetic logistic regression of substitution OCCURRENCE.
#
#   Companion to figure 1D. 1D models log(dN) and log(dS), which are undefined
#   on a branch carrying no change of that type, so it necessarily excludes
#   those branches — 61% of them for the whole gene and 99% for the ridge.
#   Reviewer 1 objected to that exclusion. This figure keeps every branch by
#   asking a question a zero can answer: did any such change occur at all?
#
#   Together the two are a hurdle: supp9 models WHETHER change happened, 1D
#   models HOW FAST it went given that it did.
#
#   Source: 04_epitope_analysis/17_phyloglm_slac_occurrence.R
#     any_N ~ 6 climate + log(EN x t) + as.factor(year)
#     any_S ~ 6 climate + log(ES x t) + as.factor(year)
#
#   Coefficients are on the log-odds scale: exp(b) is the odds ratio per one
#   standard deviation of the predictor (about 10 C for temperature). The
#   exposure term log(sites x branch duration) is in the model but not plotted —
#   it accounts for the fact that a longer branch with more scoreable codons is
#   mechanically more likely to carry a change, and is a measurement construct
#   rather than a hypothesis.
#
#   Delta AIC compares the full model against one keeping exposure and the year
#   factor but dropping the climate block, so it reads as "what do the climate
#   predictors add".
#
# INPUT:  revisions/data/rates/phylolm/all_phyloglm_occurrence_results.csv
# OUTPUT: revisions/paper/figures/supp9_phyloglm_occurrence.png / .pdf

suppressPackageStartupMessages({
  library(here)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
})

source(here("revisions", "config.R"))
source(here("revisions", "R", "forest_utils.R"))

out_dir <- here("revisions", "paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS  <- 9
BS2 <- 10

PRED <- names(PREDICTOR_LABELS)

RESP_LABELS <- c(any_N = "Non-synonymous change occurred",
                 any_S = "Synonymous change occurred")

res_path <- file.path(PATHS$rates_phylolm, "all_phyloglm_occurrence_results.csv")
if (!file.exists(res_path))
  stop("Not found: ", basename(res_path),
       "\nRun 04_epitope_analysis/17_phyloglm_slac_occurrence.R first.")

raw <- read_csv(res_path, show_col_types = FALSE) |>
  filter(term %in% PRED, gene %in% c("H1", "H3"),
         response %in% names(RESP_LABELS)) |>
  mutate(
    subtype   = ifelse(gene == "H1", "H1N1", "H3N2"),
    region    = factor(REGION_LABELS[subset], levels = REGION_LEVELS),
    resp_lab  = factor(RESP_LABELS[response], levels = RESP_LABELS),
    predictor = factor(PREDICTOR_LABELS[term], levels = PREDICTOR_LABELS),
    ci_lo     = estimate - 1.96 * std_error,
    ci_hi     = estimate + 1.96 * std_error,
    sig_fdr   = !is.na(p.fdr) & p.fdr < 0.05
  ) |>
  filter(!is.na(region))

if (nrow(raw) == 0) stop("No plottable rows in ", basename(res_path))

raw$point_fill <- ifelse(raw$sig_fdr,
                         PREDICTOR_COLOURS[as.character(raw$predictor)], "white")

# Delta AIC per panel cell; absent if script 17 ran with --skip-null.
daic <- raw |>
  distinct(resp_lab, subtype, region, AIC, AIC_null) |>
  filter(!is.na(AIC_null)) |>
  mutate(txt = sprintf("ΔAIC = %.1f", AIC_null - AIC))

p <- ggplot(raw, aes(x = estimate, y = region,
                     colour = predictor, group = predictor)) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey50",
             linewidth = 0.4) +
  geom_linerange(aes(xmin = ci_lo, xmax = ci_hi), linewidth = 0.6,
                 position = position_dodge(width = 0.7)) +
  geom_point(aes(fill = point_fill), shape = 21, size = 2.1, stroke = 0.7,
             position = position_dodge(width = 0.7))

if (nrow(daic) > 0)
  p <- p + geom_text(data = daic, aes(x = Inf, y = region, label = txt),
                     hjust = 1.05, inherit.aes = FALSE, size = 2.2,
                     colour = "grey40")

p <- p +
  scale_colour_manual(values = PREDICTOR_COLOURS, name = NULL) +
  scale_fill_identity() +
  scale_y_discrete(limits = rev, drop = TRUE) +
  scale_x_continuous(expand = expansion(mult = c(0.06, 0.34))) +
  facet_grid(resp_lab ~ subtype, scales = "free_x") +
  labs(
    title    = "Occurrence of substitution versus climate, with zeroes retained",
    subtitle = paste0("Phylogenetic logistic regression (phyloglm). Log-odds ",
                      "coefficients ± 95% CI; filled = BH-FDR p < 0.05.\n",
                      "Exposure log(sites × branch duration) is fitted but ",
                      "not shown. ΔAIC is against a null keeping exposure ",
                      "and year."),
    x = "Log-odds coefficient", y = NULL
  ) +
  theme_bw(base_size = BS) +
  theme(
    plot.title       = element_text(size = BS2, face = "bold"),
    plot.subtitle    = element_text(size = BS - 1.5, colour = "grey30"),
    strip.background = element_blank(),
    strip.text       = element_text(face = "bold", size = BS - 0.5),
    legend.position  = "bottom",
    legend.text      = element_text(size = BS - 1),
    legend.key.size  = unit(0.6, "lines"),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  ) +
  guides(colour = guide_legend(nrow = 2))

ggsave(file.path(out_dir, "supp9_phyloglm_occurrence.png"), p,
       width = 8.27, height = 6.2, dpi = 300)
ggsave(file.path(out_dir, "supp9_phyloglm_occurrence.pdf"), p,
       width = 8.27, height = 6.2, device = grDevices::cairo_pdf)

message("Written: supp9_phyloglm_occurrence.png / .pdf")

message("\nTemperature odds ratios (exp(b) per ~10 C):")
raw |>
  filter(term == "scaled_temp") |>
  mutate(OR = round(exp(estimate), 3), estimate = round(estimate, 4),
         p.fdr = signif(p.fdr, 3),
         sig = ifelse(sig_fdr, "*", "")) |>
  select(gene, region, response, n, n_events, estimate, OR, p.fdr, sig) |>
  arrange(response, gene, region) |>
  as.data.frame() |>
  print(row.names = FALSE)
