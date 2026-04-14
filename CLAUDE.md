# CLAUDE.md — Context for AI assistants

## Project purpose

This repository analyses how environmental variables (temperature, relative humidity) influence the molecular evolution of Influenza A surface proteins. Four gene segments are studied: **H1, H3, N1, N2**.

The pipeline takes manually exported GISAID sequences + metadata, attaches ERA5 climate data, builds phylogenetic trees, estimates dN/dS via PAML, dates the tree with TreeTime, and tests environmental associations using phylogenetic comparative regression (phylolm).

---

## Data

- **GISAID files** live in `data/raw/` and are **gitignored**. They must be manually exported from gisaid.org. Never attempt to download them programmatically.
- Expected filenames: `data/raw/{gene}_sequences.fasta` and `data/raw/{gene}_metadata.tsv`
- All processed data (metadata CSVs, alignments, trees, rates) are also gitignored — only code is tracked.

---

## Code conventions

### Paths — never hardcode
All paths use `here::here()` via `config.R`. Every R script begins with `source(here("config.R"))`. The central path list is `PATHS`, tool paths are `TOOLS`, and analysis parameters are `PARAMS`.

### One script per stage, gene as parameter
Scripts loop over `GENES <- c("H1", "H3", "N1", "N2")` from `config.R`. Each script accepts an optional `--gene X` command-line argument to run a single segment. Do not create per-gene copies of scripts.

### R functions in `R/`
All reusable logic lives in `R/*.R` utility files with roxygen2-style comments. Scripts in `scripts/` should be thin wrappers that call these functions.

### Error handling
External tool calls (`system2()`) should check return codes or output files. Failed pairwise batches (MACSE, codeml) are logged to `failed_pairs.txt` in the output directory and processing continues — do not stop the whole pipeline on a single failure.

### Python scripts
Use `argparse` for all CLI arguments. Check that GEE is authenticated before processing. ERA5 temperature values are in Kelvin — always convert to Celsius.

---

## External tools

Paths are set in `config.R` under `TOOLS`. Defaults use `Sys.which()` where possible. MACSE and MAPLE require explicit paths.

| Variable | Tool |
|----------|------|
| `TOOLS$mafft` | MAFFT (multiple alignment) |
| `TOOLS$macse` | MACSE v2.07 JAR path |
| `TOOLS$maple` | MAPLE executable |
| `TOOLS$codeml` | PAML codeml |
| `TOOLS$treetime` | TreeTime |

---

## Key scientific decisions

- **phylolm, not brms**: Frequentist phylogenetic regression is sufficient and faster; Bayesian hierarchical models are not used in this version.
- **ERA5-Land**: Chosen for spatial resolution (~0.1°) and temporal coverage (daily aggregates back to 1950). Monthly means are fetched for each sample's lat/lon/date.
- **Pairwise MACSE**: Codon alignment is done pairwise (batch × batch) rather than all-vs-all to keep runtime manageable. Batches are labelled A, B, C... and all pairwise combinations (AA, AB, ...) are processed.
- **Two-round MAPLE pruning**: Round 1 removes very short branches (likely duplicate/identical sequences), round 2 removes very long branches (likely sequencing errors or mis-assignments). Thresholds in `PARAMS`.
- **Relative humidity**: Derived from ERA5 temperature + dewpoint using the Magnus formula (see `01_attach_climate.py`).

---

## Stage 00 — metadata processing

The metadata preparation scripts (`scripts/00_metadata/`) are **subject to change** as methods improve. Keep them modular and well-commented. When updating, ensure:

1. Output column names stay consistent (`accession`, `latitude`, `longitude`, `collection_date`, `year`, `month`, `country`, `region`, `host`) so downstream scripts don't break.
2. The climate join key is always `accession`.

---

## TODOs / known limitations

- Population density data is not yet integrated (columns `pop_density` / `scaled_pop_density` will be NA until a source is added).
- Urban/rural subset thresholds in `02_phylolm.R` are provisional (> 500 people/km² = urban).
- MACSE pairwise processing is the computational bottleneck — consider parallelisation with `parallel::mclapply` or a job scheduler for large datasets.
