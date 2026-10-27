# scripts/04_epitope_analysis/08_terminal_vs_internal.R
#
# PURPOSE: Sanity check on how much of the terminal-branch signal is weakly
#          deleterious. Compares pooled omega on TERMINAL branches against
#          INTERNAL branches, per gene and region.
#
# THE ARGUMENT: a non-synonymous mutation that is weakly deleterious can be
#   observed on a terminal branch — it arose recently and purifying selection
#   has not yet removed it — but will rarely survive to sit on an internal
#   branch, which by definition leads to descendants that persisted. So
#
#     terminal omega >> internal omega   suggests a large transient deleterious
#                                        load in the terminal-branch data
#     terminal omega ~  internal omega   suggests terminal changes are broadly
#                                        the same material that goes on to fix
#
#   This matters because figures 1D and 2 are built entirely on terminal
#   branches. The ratio per region says where that load is concentrated.
#
# WHY POOLED, NOT AVERAGED: per-branch omega is degenerate here — with N and S
#   usually 0 or 1 on a terminal branch, omega collapses to (N/S) times a
#   near-constant, and is undefined whenever S = 0 (97%+ of ridge branches).
#   Averaging those would be meaningless. Instead counts and expected sites are
#   summed across branches first and the ratio taken once:
#
#     omega = (sum N / sum EN) / (sum S / sum ES)
#
#   which is the standard aggregate estimator and uses every branch, including
#   those carrying no change at all.
#
# WHY OMEGA AND NOT dN: internal and terminal branches differ in length, so dN
#   and dS are not comparable between them without normalising by branch
#   duration. Omega is a ratio over the same branch, so length cancels and the
#   comparison is valid as it stands. dN and dS are still written to the CSV for
#   reference, but they should NOT be compared across branch types.
#
# UNCERTAINTY: percentile bootstrap over branches (--boot replicates, default
#   1000), resampling branches within each gene x region x branch type. This
#   matters most for the ridge, where the whole internal-branch estimate rests
#   on a few dozen non-synonymous events.
#
# CAVEATS to carry into the discussion:
#   - Internal branches are not uniformly long in a densely sampled tree; this
#     is not a clean long-vs-short contrast.
#   - Longer internal branches carry slightly more multiple-hit saturation,
#     which depresses internal omega and so inflates the ratio. Small at these
#     divergences, but directional, making the ratios mild upper bounds.
#   - Internal branches are shared ancestry, so they are not independent of the
#     terminal branches descending from them.
#
# INPUT:
#   data/rates/{gene}_slac_branches.csv   from 03_slac_dnds.R
#   data/slac/{gene}/tree.nwk             to identify terminal nodes
#
# OUTPUT:
#   data/rates/terminal_vs_internal_dnds.csv
#   figures/epitope/terminal_vs_internal_omega.png / .pdf
#
# FLAGS:
#   --gene X     run a single gene (default: H1 H3)
#   --boot N     bootstrap replicates (default 1000; 0 to skip)
#   --seed N     default 42

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
  library(viridisLite)
})

source(here("config.R"))
source(here("R", "forest_utils.R"))

args  <- commandArgs(trailingOnly = TRUE)
genes <- if ("--gene" %in% args) args[which(args == "--gene") + 1] else c("H1", "H3")
n_boot <- if ("--boot" %in% args) as.integer(args[which(args == "--boot") + 1]) else 1000L
seed   <- if ("--seed" %in% args) as.integer(args[which(args == "--seed") + 1]) else 42L

set.seed(seed)

REGIONS <- c("full", "HA1", "HA2", "epitope", "ridge")
TYPE_COLOURS <- setNames(viridis(2, end = 0.75), c("Internal", "Terminal"))

fig_dir <- here("figures", "epitope")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


#' Pooled omega from vectors of per-branch counts and expected sites.
pooled_omega <- function(N, S, EN, ES, idx = NULL) {
  if (!is.null(idx)) { N <- N[idx]; S <- S[idx]; EN <- EN[idx]; ES <- ES[idx] }
  sEN <- sum(EN, na.rm = TRUE); sES <- sum(ES, na.rm = TRUE)
  if (sEN <= 0 || sES <= 0) return(NA_real_)
  dN <- sum(N, na.rm = TRUE) / sEN
  dS <- sum(S, na.rm = TRUE) / sES
  if (!is.finite(dS) || dS <= 0) return(NA_real_)
  dN / dS
}

#' Percentile bootstrap CI for pooled omega, resampling branches.
boot_omega <- function(N, S, EN, ES, reps) {
  if (reps < 1 || length(N) < 2) return(c(NA_real_, NA_real_))
  n <- length(N)
  vals <- vapply(seq_len(reps), function(i)
    pooled_omega(N, S, EN, ES, sample.int(n, n, replace = TRUE)), numeric(1))
  unname(stats::quantile(vals, c(0.025, 0.975), na.rm = TRUE))
}

# ── Compute ───────────────────────────────────────────────────────────────────

rows <- list()

for (gene in genes) {
  br_path   <- file.path(PATHS$rates, paste0(gene, "_slac_branches.csv"))
  tree_path <- file.path(here(), "data", "slac", gene, "tree.nwk")
  if (!file.exists(br_path)) {
    warning(gene, ": ", basename(br_path), " not found — run 03_slac_dnds.R"); next
  }
  if (!file.exists(tree_path)) { warning(gene, ": ", tree_path, " not found"); next }

  message("\n[", gene, "]")
  b  <- suppressWarnings(read_csv(br_path, show_col_types = FALSE,
                                  guess_max = 100000))
  tr <- read.tree(tree_path)
  b$type <- ifelse(b$node %in% tr$tip.label, "Terminal", "Internal")

  for (r in intersect(REGIONS, unique(b$region))) {
    for (ty in c("Internal", "Terminal")) {
      d <- b[b$region == r & b$type == ty, , drop = FALSE]
      d <- d[!is.na(d$EN) & !is.na(d$ES), , drop = FALSE]
      if (nrow(d) < 2) next

      om <- pooled_omega(d$N, d$S, d$EN, d$ES)
      ci <- boot_omega(d$N, d$S, d$EN, d$ES, n_boot)

      rows[[paste(gene, r, ty)]] <- tibble(
        gene = gene, region = r, type = ty,
        branches = nrow(d),
        N = sum(d$N, na.rm = TRUE), S = sum(d$S, na.rm = TRUE),
        EN = sum(d$EN, na.rm = TRUE), ES = sum(d$ES, na.rm = TRUE),
        dN = sum(d$N, na.rm = TRUE) / sum(d$EN, na.rm = TRUE),
        dS = sum(d$S, na.rm = TRUE) / sum(d$ES, na.rm = TRUE),
        omega = om, omega_lo = ci[1], omega_hi = ci[2]
      )
    }
    oi <- rows[[paste(gene, r, "Internal")]]$omega
    ot <- rows[[paste(gene, r, "Terminal")]]$omega
    if (!is.null(oi) && !is.null(ot))
      message(sprintf("  %-8s internal %.4f   terminal %.4f   ratio %.2f",
                      r, oi, ot, ot / oi))
  }
}

if (length(rows) == 0) stop("No results computed.")
res <- bind_rows(rows)

# Ratio and its bootstrap-free interval (CIs are on each omega separately, so
# the ratio interval below is indicative rather than exact).
ratio <- res |>
  select(gene, region, type, omega, omega_lo, omega_hi) |>
  pivot_wider(names_from = type, values_from = c(omega, omega_lo, omega_hi)) |>
  mutate(ratio = omega_Terminal / omega_Internal,
         ratio_lo = omega_lo_Terminal / omega_hi_Internal,
         ratio_hi = omega_hi_Terminal / omega_lo_Internal)

res <- res |> left_join(ratio |> select(gene, region, ratio, ratio_lo, ratio_hi),
                        by = c("gene", "region"))

out_csv <- file.path(PATHS$rates, "terminal_vs_internal_dnds.csv")
write_csv(res, out_csv)
message("\nWritten: ", basename(out_csv))

message("\nPooled omega, terminal vs internal:")
ratio |>
  mutate(across(c(omega_Internal, omega_Terminal), \(x) round(x, 4)),
         ratio = round(ratio, 2),
         ratio_CI = sprintf("[%.2f, %.2f]", ratio_lo, ratio_hi)) |>
  select(gene, region, omega_Internal, omega_Terminal, ratio, ratio_CI) |>
  as.data.frame() |>
  print(row.names = FALSE)

# ── Figure ────────────────────────────────────────────────────────────────────

# One panel. The earlier two-panel version put the ratio in its own panel, but
# the ratio is just the gap between the two points here, so it is annotated
# instead. y is numeric rather than discrete so the two branch types can be
# offset to their own rows — on a shared row their confidence intervals overlap
# and become unreadable.
Y_OFFSET <- 0.17
y_levels <- rev(REGION_LEVELS)          # position 1 is the bottom of the panel

plot_df <- res |>
  mutate(region_lab = REGION_LABELS[region],
         subtype    = ifelse(gene %in% c("H1", "N1"), "H1N1", "H3N2"),
         type       = factor(type, levels = c("Internal", "Terminal"))) |>
  filter(region_lab %in% y_levels) |>
  mutate(y_base = match(region_lab, y_levels),
         y      = y_base + ifelse(type == "Terminal", Y_OFFSET, -Y_OFFSET))

# Connector showing the gap between the two estimates within a region.
seg <- plot_df |>
  select(subtype, y_base, type, omega) |>
  pivot_wider(names_from = type, values_from = omega) |>
  mutate(y_int = y_base - Y_OFFSET, y_term = y_base + Y_OFFSET)

# Ratio label, placed at the right edge. Dark where the bootstrap CI for the
# ratio excludes 1, grey where it does not — so "no evidence of excess" reads
# differently from "excess".
rat_df <- ratio |>
  mutate(region_lab = REGION_LABELS[region],
         subtype    = ifelse(gene %in% c("H1", "N1"), "H1N1", "H3N2")) |>
  filter(region_lab %in% y_levels) |>
  mutate(y_base   = match(region_lab, y_levels),
         excludes = !is.na(ratio_lo) & !is.na(ratio_hi) &
                    (ratio_lo > 1 | ratio_hi < 1),
         txt      = sprintf("%.2f×", ratio),
         col      = ifelse(excludes, "grey15", "grey65"))

p <- ggplot(plot_df, aes(x = omega, y = y)) +
  geom_segment(data = seg,
               aes(x = Internal, xend = Terminal, y = y_int, yend = y_term),
               inherit.aes = FALSE, colour = "grey75", linewidth = 0.7) +
  geom_linerange(aes(xmin = omega_lo, xmax = omega_hi, colour = type),
                 linewidth = 0.8, show.legend = FALSE) +
  geom_point(aes(colour = type), size = 2.8) +
  # show.legend = FALSE: without it this layer contributes its own glyph ("a")
  # to the colour guide alongside Internal/Terminal.
  geom_text(data = rat_df,
            aes(x = Inf, y = y_base, label = txt, colour = I(col)),
            hjust = 1.15, inherit.aes = FALSE, size = 3,
            show.legend = FALSE) +
  scale_colour_manual(values = TYPE_COLOURS, name = NULL) +
  scale_y_continuous(breaks = seq_along(y_levels), labels = y_levels,
                     limits = c(0.4, length(y_levels) + 0.6)) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.20))) +
  facet_wrap(~ subtype, scales = "free_x") +
  labs(title = "Pooled dN/dS on internal vs terminal branches",
       subtitle = paste0(
         "Aggregate estimator (ΣN/ΣEN)/(ΣS/ΣES); bars are ",
         n_boot, "-replicate bootstrap 95% CI.\n",
         "Labels give terminal/internal ω; dark where the ratio CI ",
         "excludes 1. Values above 1 mean non-synonymous\nchanges are visible ",
         "at the tips that do not persist on internal branches — a ",
         "transient weakly deleterious load."),
       x = "ω  (dN/dS)", y = NULL) +
  theme_bw(base_size = 11) +
  theme(legend.position    = "bottom",
        plot.subtitle      = element_text(size = 8.5, colour = "grey30"),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold"),
        panel.grid.minor   = element_blank(),
        panel.grid.major.y = element_blank())

ggsave(file.path(fig_dir, "terminal_vs_internal_omega.png"), p,
       width = 9, height = 5.2, dpi = 300)
ggsave(file.path(fig_dir, "terminal_vs_internal_omega.pdf"), p,
       width = 9, height = 5.2, device = grDevices::cairo_pdf)
message("Saved: terminal_vs_internal_omega.png / .pdf")

message("\nDone.")
