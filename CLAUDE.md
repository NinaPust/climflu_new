# CLAUDE.md — Context for AI assistants

## Project purpose

This repository analyses how environmental variables (temperature, relative humidity) influence the molecular evolution of Influenza A surface proteins. Six gene segments are studied: **H1, H3, N1, N2** (surface proteins) plus **PA_H1N1, PA_H3N2** (internal protein, negative control).

The pipeline takes manually exported GISAID sequences + metadata, attaches ERA5 climate data, builds phylogenetic trees, estimates dN/dS via PAML, dates the tree with TreeTime, and tests environmental associations using phylogenetic comparative regression (phylolm).

---

## Data

- **GISAID files** live in `data/raw/` and are **gitignored**. They must be manually exported from gisaid.org. Never attempt to download them programmatically.
- Raw FASTA files may be split by date range (e.g. `H3_2000_18.fasta` + `H3_2019.fasta`). The parse script merges these into canonical `{gene}_sequences.fasta` files.
- Raw metadata is in XLS format (GISAID EpiFlu export), not TSV. One file per subtype: `H1N1_meta.xls`, `H3N2_2000_18_meta.xls` + `H3N2_2019_metadata.xls`.
- GISAID EpiFlu XLS exports do **not** include lat/lon. Coordinates are added by geocoding the Location string via Nominatim (script `01_geocode.R`), with results cached in `data/raw/location_coords_cache.csv`.
- All processed data (metadata CSVs, alignments, trees, rates) are also gitignored — only code is tracked.

---

## Code conventions

### Paths — never hardcode
All paths use `here::here()` via `config.R`. Every R script begins with `source(here("config.R"))`. The central path list is `PATHS`, tool paths are `TOOLS`, and analysis parameters are `PARAMS`.

### One script per stage, gene as parameter
Scripts loop over `GENES <- c("H1", "H3", "N1", "N2")` from `config.R`. PA segments (`PA_H1N1`, `PA_H3N2`) are defined as `ALL_GENES` in scripts that handle them. Each script accepts an optional `--gene X` command-line argument to run a single segment. Do not create per-gene copies of scripts.

### R functions in `R/`
All reusable logic lives in `R/*.R` utility files with roxygen2-style comments. Scripts in `scripts/` should be thin wrappers that call these functions.

### Error handling
External tool calls (`system2()`) should check return codes or output files. Failed pairwise batches (MACSE, codeml) are logged to `failed_pairs.txt` in the output directory and processing continues — do not stop the whole pipeline on a single failure.

### Python scripts
Use `argparse` for all CLI arguments. ERA5 temperature values are in Kelvin — always convert to Celsius. Python scripts run in the `climflu` conda environment (miniforge3).

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
- **ERA5-Land**: Chosen for spatial resolution (~0.1°) and temporal coverage (back to 1950). Monthly means are fetched from the Copernicus Climate Data Store (CDS) via `cdsapi`. One NetCDF per year-month, point values extracted with `xarray`.
- **Pairwise MACSE**: Codon alignment is done pairwise (batch × batch) rather than all-vs-all to keep runtime manageable. Batches are labelled A, B, C... and all pairwise combinations (AA, AB, ...) are processed.
- **Two-round MAPLE pruning**: Round 1 removes very short branches (likely duplicate/identical sequences), round 2 removes very long branches (likely sequencing errors or mis-assignments). Thresholds in `PARAMS`.
- **Relative humidity**: Derived from ERA5 temperature + dewpoint using the Magnus formula (see `02_attach_climate.py`).
- **PA as negative control**: PA (polymerase acidic protein) evolves slowly and should show no climate association — serves as a methodological control for H1/H3/N1/N2.
- **Geocoding**: Location strings geocoded via Nominatim (OpenStreetMap) with fallback to coarser geographic levels. Results cached to avoid re-querying.

---

## Stage 00 — metadata processing

The metadata preparation scripts (`scripts/00_metadata/`) are **subject to change** as methods improve. Current pipeline:

1. `00_parse_gisaid.R` — reads XLS metadata + merges split FASTAs → `data/metadata/{gene}_metadata.csv`
2. `01_geocode.R` — geocodes unique location strings via Nominatim → fills lat/lon in metadata CSVs
3. `02_attach_climate.py` — fetches ERA5-Land monthly means from CDS → `data/climate/{gene}_climate.csv`

Exploratory diagnostics (gitignored): `00z_explore_metadata.R` → `figures/metadata_diagnostics/`

When updating, ensure:

1. Output column names stay consistent (`accession`, `latitude`, `longitude`, `collection_date`, `year`, `month`, `country`, `region`, `host`) so downstream scripts don't break.
2. The climate join key is `(latitude, longitude, year, month)`.

---

## TODOs / known limitations

- Population density data is not yet integrated (columns `pop_density` / `scaled_pop_density` will be NA until a source is added).
- Urban/rural subset thresholds in `02_phylolm.R` are provisional (> 500 people/km² = urban).
- MACSE pairwise processing is the computational bottleneck — consider parallelisation with `parallel::mclapply` or a job scheduler for large datasets.
