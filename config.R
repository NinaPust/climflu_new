# config.R
# Central configuration for the climflu_new analysis pipeline.
# All scripts source this file. Do NOT hardcode paths in scripts.

library(here)

# ── Gene segments ─────────────────────────────────────────────────────────────
GENES <- c("H1", "H3", "N1", "N2")

# ── Data paths (relative to project root via here::here()) ───────────────────
PATHS <- list(
  raw        = here("data", "raw"),
  metadata   = here("data", "metadata"),
  climate    = here("data", "climate"),
  population = here("data", "population"),
  alignments = here("data", "alignments"),
  trees      = here("data", "trees"),
  trees_maple  = here("data", "trees", "raw", "maple"),
  trees_fast   = here("data", "trees", "raw", "fasttree"),
  trees_pruned = here("data", "trees", "pruned"),
  trees_timed  = here("data", "trees", "timed"),

  # Test subset paths (mirror of main structure, populated by 00_subsample.R)
  test_trees_pruned = here("data", "test", "trees", "pruned"),
  test_trees_timed  = here("data", "test", "trees", "timed"),
  test_metadata     = here("data", "test", "metadata"),
  test_climate      = here("data", "test", "climate"),
  test_rates        = here("data", "test", "rates"),
  codeml     = here("data", "codeml"),
  yn00       = here("data", "yn00"),
  rates      = here("data", "rates"),
  scripts    = here("scripts"),
  R          = here("R")
)

# ── External tool paths ───────────────────────────────────────────────────────
# Update these to match the tool locations on your system.
TOOLS <- list(
  mafft    = Sys.which("mafft"),       # or full path, e.g. "/usr/local/bin/mafft"
  macse    = "/path/to/macse_v2.07.jar",
  maple    = "/home/nina/influenza_env/programs/maple",
  codeml   = "/home/nina/paml4.8/bin/codeml",
  yn00     = "/home/nina/paml4.8/bin/yn00",
  treetime = "/home/nina/miniforge3/envs/treetime_env/bin/treetime"
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

message("config.R loaded. Project root: ", here())
