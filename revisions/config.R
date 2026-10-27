# revisions/config.R
# Isolated configuration for the manuscript-revision pipeline.
# Mirrors the project-root config.R, but all PATHS resolve under revisions/,
# so revision work never reads or writes the original submitted analysis.
#
# Paths for pipeline stages unaffected by the reviewer-requested fixes (raw,
# metadata, climate, population, alignments, trees, codeml, yn00, test) resolve
# through symlinks in revisions/data/ back to the original data/ directories —
# read-only reuse, no duplication. Paths for stages that get re-run with fixes
# (epitope, rates) point to fresh directories that live only under revisions/.
#
# All scripts under revisions/scripts/ should source(here("revisions", "config.R"))
# instead of source(here("config.R")).

library(here)

# ── Gene segments ─────────────────────────────────────────────────────────────
GENES <- c("H1", "H3", "N1", "N2")

# ── Data paths (relative to project root via here::here()) ───────────────────
PATHS <- list(
  raw        = here("revisions", "data", "raw"),
  metadata   = here("revisions", "data", "metadata"),
  climate    = here("revisions", "data", "climate"),
  population = here("revisions", "data", "population"),
  alignments = here("revisions", "data", "alignments"),
  trees      = here("revisions", "data", "trees"),
  trees_maple  = here("revisions", "data", "trees", "raw", "maple"),
  trees_fast   = here("revisions", "data", "trees", "raw", "fasttree"),
  trees_pruned = here("revisions", "data", "trees", "pruned"),
  trees_timed  = here("revisions", "data", "trees", "timed"),

  # Test subset paths (mirror of main structure, populated by 00_subsample.R)
  test_trees_pruned = here("revisions", "data", "test", "trees", "pruned"),
  test_trees_timed  = here("revisions", "data", "test", "trees", "timed"),
  test_metadata     = here("revisions", "data", "test", "metadata"),
  test_climate      = here("revisions", "data", "test", "climate"),
  test_rates        = here("revisions", "data", "test", "rates"),
  codeml     = here("revisions", "data", "codeml"),
  yn00       = here("revisions", "data", "yn00"),
  epitope    = here("revisions", "data", "epitope"),
  rates         = here("revisions", "data", "rates"),
  rates_phylolm = here("revisions", "data", "rates", "phylolm"),
  rates_brms    = here("revisions", "data", "rates", "brms"),
  scripts    = here("revisions", "scripts"),
  R          = here("revisions", "R")
)

# ── External tool paths ───────────────────────────────────────────────────────
# Update these to match the tool locations on your system.
TOOLS <- list(
  mafft    = Sys.which("mafft"),       # or full path, e.g. "/usr/local/bin/mafft"
  macse    = "/path/to/macse_v2.07.jar",
  maple    = "/home/nina/influenza_env/programs/maple",
  codeml   = "/home/nina/paml4.8/bin/codeml",
  yn00     = "/home/nina/paml4.8/bin/yn00",
  treetime = "/home/nina/miniforge3/envs/treetime_env/bin/treetime",
  hyphy    = "/home/nina/miniforge3/envs/hyphy_env/bin/hyphy"
)

# Warn if any tool is not found
missing_tools <- names(TOOLS)[sapply(TOOLS, function(p) nchar(p) == 0 || p == "/path/to/MAPLE" || grepl("^/path/to/", p))]
if (length(missing_tools) > 0) {
  warning("Tool paths not configured in config.R: ", paste(missing_tools, collapse = ", "))
}

# ── Analysis parameters ───────────────────────────────────────────────────────
PARAMS <- list(
  # Date range for sample inclusion
  date_from = "2009-01-01",
  date_to   = "2023-12-31",

  # MAPLE tree pruning
  prune_min_bl = 0.0, # Fixed lower cutoff: removes near-identical duplicates
  prune_long_k = 4,   # MAD multiplier for long-branch outlier detection

  # MACSE batch size (number of accessions per lettered subset)
  macse_batch_size = 50,

  # Minimum gap threshold for alignment trimming (codon columns)
  trim_gap_threshold = 0.5,

  # PAML codeml model (0 = one ratio, 1 = free ratios)
  codeml_model      = 0,
  codeml_batch_size = 100,  # sequences per batch; ~80 batches for 8k seqs → ~3,240 pairs
  codeml_cores      = 16,   # parallel jobs; ~25% of 64-core server (polite default)
  codeml_nice       = 10,   # nice priority for codeml processes (0=normal, 19=lowest)

  # phylolm: Pagel's lambda method
  phylolm_method = "lambda",

  # Number of threads for MAFFT
  mafft_threads = 4
)

message("revisions/config.R loaded. Project root: ", here())
