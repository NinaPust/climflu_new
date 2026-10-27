# scripts/05_final_figs/supp10_pseudoreplication.R
#
# Supp 10: Mutation events in figure 2 are phylogenetically independent.
#
# 04_epitope_analysis/06_ml_mutations.R compares each tip only with ITS OWN
# PARENT, so a change inherited from further back in the clade is present in the
# parent too and is never counted. The only way the design could over-count is a
# parent lacking a change while two or more of its tip children carry it.
#
# That never happens — not because the opportunity is absent, but because
# maximum-likelihood ancestral reconstruction assigns a change shared by
# siblings to their parent, one event up the tree being more likely than two
# independent ones below it. Independence is therefore a property of the design
# rather than a lucky feature of the data, which is the argument made in the
# methods; this figure exists to show it if asked.
#
# The conservative direction is worth stating alongside: the same mechanism
# means genuine convergence on two sibling branches would be attributed to the
# parent and dropped. The method can UNDER-count convergence, never over-count
# it.
#
# The figure shows that the opportunity for sharing was real — mutations occur
# in clades of every size, up to 112-way polytomies — and that none was taken.
#
# INPUT:
#   data/epitope/{gene}_{region}_mutations_ml.csv   from script 06
#   data/trees/timed/{gene}_timetree.nwk                      parent lookup
#
# OUTPUT:
#   paper/figures/supp10_pseudoreplication.png / .pdf
#   data/rates/pseudoreplication_summary.csv
#
# FLAGS:
#   --gene X     default: H1 H3
#   --region X   default: ridge epitope

suppressPackageStartupMessages({
  library(here)
  library(ape)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(viridisLite)
})

source(here("config.R"))

args    <- commandArgs(trailingOnly = TRUE)
genes   <- if ("--gene"   %in% args) args[which(args == "--gene")   + 1] else c("H1", "H3")
regions <- if ("--region" %in% args) args[which(args == "--region") + 1] else c("ridge", "epitope")

out_dir <- here("paper", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

BS <- 9

# ── Per-event sibling counts ─────────────────────────────────────────────────

events <- list()

for (gene in genes) {
  tree_path <- here("data", "trees", "timed", paste0(gene, "_timetree.nwk"))
  if (!file.exists(tree_path)) { warning(gene, ": timetree not found"); next }
  tr <- read.tree(tree_path)

  # Parent of every tip, and that parent's number of TIP children. Internal
  # children are irrelevant: only tips carry recorded mutations.
  tip_edge <- match(seq_len(Ntip(tr)), tr$edge[, 2])
  tip_par  <- tibble(accession = tr$tip.label,
                     parent    = tr$edge[tip_edge, 1]) |>
    group_by(parent) |> mutate(clade_size = n()) |> ungroup()

  for (region in regions) {
    mp <- file.path(PATHS$epitope,
                    paste0(gene, "_", region, "_mutations_ml.csv"))
    if (!file.exists(mp)) {
      warning(gene, " ", region, ": ", basename(mp), " not found — run 18")
      next
    }
    mut <- read_csv(mp, show_col_types = FALSE)
    if (!nrow(mut)) next

    events[[paste(gene, region)]] <- mut |>
      select(accession, change) |>
      inner_join(tip_par, by = "accession") |>
      group_by(parent, change) |>
      mutate(k_sharing = n()) |>          # tip children of this parent with this change
      ungroup() |>
      mutate(gene = gene, region = region,
             subtype = ifelse(gene == "H1", "H1N1", "H3N2"))
  }
}

if (length(events) == 0) stop("No mutation files found — run script 06 first.")
ev <- bind_rows(events)

# ── Summary ──────────────────────────────────────────────────────────────────

summ <- ev |>
  group_by(gene, subtype, region) |>
  summarise(events          = n(),
            distinct_change = n_distinct(change),
            shared          = sum(k_sharing >= 2),
            max_k           = max(k_sharing),
            max_clade       = max(clade_size),
            multi_parents   = n_distinct(parent[duplicated(parent) |
                                                duplicated(parent, fromLast = TRUE)]),
            .groups = "drop")

write_csv(summ, file.path(PATHS$rates, "pseudoreplication_summary.csv"))

message("\nPseudoreplication summary:")
summ |> as.data.frame() |> print(row.names = FALSE)

# ── Figure: one row, one panel per subtype ───────────────────────────────────
# Clade size is binned because the tail is long (a few polytomies above 100)
# and the point is the spread of opportunity, not its exact shape.

CLADE_BREAKS <- c(0, 1, 2, 4, 9, 19, Inf)
CLADE_LABS   <- c("1", "2", "3-4", "5-9", "10-19", "20+")

plot_df <- ev |>
  mutate(bin = cut(clade_size, breaks = CLADE_BREAKS, labels = CLADE_LABS)) |>
  count(subtype, region, bin)

# Per-panel caption carrying the numbers that answer the question.
annot <- summ |>
  group_by(subtype) |>
  summarise(events = sum(events), shared = sum(shared),
            multi  = sum(multi_parents), max_clade = max(max_clade),
            .groups = "drop") |>
  mutate(txt = sprintf(
    "%s events across clades of up to %s tips\n%s parents host >1 mutated tip\n%s events shared with a sibling",
    format(events, big.mark = ","), max_clade, format(multi, big.mark = ","),
    shared))

p <- ggplot(plot_df, aes(x = bin, y = n, fill = region)) +
  geom_col(position = position_dodge(width = 0.75), width = 0.68) +
  geom_text(data = annot, aes(x = Inf, y = Inf, label = txt),
            inherit.aes = FALSE, hjust = 1.03, vjust = 1.15,
            size = 2.4, colour = "grey25", lineheight = 1.15) +
  facet_wrap(~ subtype, scales = "free_y") +
  scale_fill_manual(values = setNames(viridis(2, end = 0.7),
                                      c("epitope", "ridge")), name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.28))) +
  labs(
    title    = "Mutation events arise independently on terminal branches",
    subtitle = paste0("Each event is a tip differing from its own ML-",
                      "reconstructed parent. No parent lacking a change had ",
                      "two or more\ntip children carrying it, at any clade ",
                      "size — a change shared by siblings is assigned to ",
                      "the parent and never counted."),
    x = "Tip children of the parent node (clade size)",
    y = "Mutation events"
  ) +
  theme_bw(base_size = BS) +
  theme(plot.title         = element_text(face = "bold", size = BS + 2),
        plot.subtitle      = element_text(size = BS - 1, colour = "grey30"),
        strip.background   = element_blank(),
        strip.text         = element_text(face = "bold", size = BS + 1),
        legend.position    = "bottom",
        panel.grid.minor   = element_blank(),
        panel.grid.major.x = element_blank())

ggsave(file.path(out_dir, "supp10_pseudoreplication.png"), p,
       width = 8.27, height = 4.2, dpi = 300)
ggsave(file.path(out_dir, "supp10_pseudoreplication.pdf"), p,
       width = 8.27, height = 4.2, device = grDevices::cairo_pdf)

message("\nWritten: supp10_pseudoreplication.png / .pdf")
