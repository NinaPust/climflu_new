# climflu_new

Studying the effects of environmental variables (temperature, relative humidity) on the molecular evolution of Influenza A surface proteins.

**Gene segments:** H1, H3, N1, N2 + PA_H1N1, PA_H3N2 (PA as negative control)  
**Environmental data:** ERA5-Land (Copernicus CDS), GHSL population density  
**Statistical approach:** phylogenetic comparative regression (phylolm)

---

## Pipeline overview

```
Stage 00 — Metadata & climate
  GISAID export (manual) → parse metadata → geocode locations (Nominatim)
  → attach ERA5-Land climate (temp, dewpoint → RH) + GHSL population density

Stage 01 — Alignment
  MAFFT multiple alignment → filter (>50% gaps, duplicates removed)

Stage 02 — Phylogenetics & dating
  Round 1: MAPLE tree → TreeTime molecular dating → clock filter (±4 IQR outliers removed)
  Round 2: MAPLE tree (filtered) → TreeTime chronogram (final)
  → per-branch substitution rates = MAPLE branch length / TreeTime duration

Stage 03 — dN/dS estimation
  Parsimony ancestral reconstruction (phangorn) → count synonymous/non-synonymous
  changes per tip branch → dN, dS, ω rates
  Subregion analysis: extract HA2, epitope, antigenic ridge → same parsimony approach

Stage 04 — Statistical analysis
  phylolm regression: log rate ~ temp + RH + ΔT + ΔRH + log(pop density)
                      + log(sampling density) + year (factor)
  → Benjamini–Hochberg FDR correction
  Subsampling stability: 1,000 year-stratified subsamples (1000 seqs / 15% / 30% / 60%)
  Recurrent mutation analysis: Wilcoxon rank-sum tests for mutations ≥5 occurrences
```

---

## Data

**GISAID sequences cannot be downloaded programmatically.** You must export them manually from [gisaid.org](https://gisaid.org).

### Required exports

From the EpiFlu database, download sequences and metadata for each subtype (H1N1, H3N2) and for PA segments. Filter to **original passage only**.

| File type | Examples |
|-----------|----------|
| Sequence FASTA | `H1_2000_19.fasta`, `H3_2000_18.fasta` + `H3_2019.fasta` |
| Metadata XLS | `H1N1_meta.xls`, `H3N2_2000_18_meta.xls` + `H3N2_2019_metadata.xls` |
| PA FASTA | `PA_H1N1_2000_19.fasta`, `PA_H3N2_2000_18.fasta` + `PA_H3N2_2019.fasta` |

FASTA files split across date ranges are automatically merged by the parse script. Lat/lon are not in GISAID EpiFlu exports — coordinates are added by geocoding (script `01_geocode.R`).

These files are gitignored and must be obtained independently.

---

## Requirements

### R packages

```r
install.packages(c(
  # Core utilities
  "here", "tidyverse",
  # Data I/O and wrangling
  "readxl", "lubridate",
  # Geocoding and spatial
  "tidygeocoder", "terra", "sf", "rnaturalearth",
  # Phylogenetics
  "ape", "phangorn", "phylolm",
  # Visualisation
  "ggplot2", "patchwork", "viridis", "forcats", "ggridges", "circlize"
))
```

Or restore the full environment:

```r
# install.packages("renv")
renv::restore()
```

### Python environment

```bash
mamba create -n climflu python=3.12 pandas cdsapi xarray netcdf4 -y
conda activate climflu
```

ERA5 climate data is fetched from the Copernicus CDS. Register at [cds.climate.copernicus.eu](https://cds.climate.copernicus.eu), accept the ERA5-Land licence, and create `~/.cdsapirc` with your API key.

### External CLI tools

| Tool | Purpose | Version tested |
|------|---------|----------------|
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Multiple sequence alignment | v7.526 |
| [MAPLE](https://github.com/NicolaDM/MAPLE) | Approximate maximum-likelihood phylogeny | v0.6.8 |
| [TreeTime](https://github.com/neherlab/treetime) | Molecular clock dating | ≥ 0.9 |

---

## Configuration

Edit `config.R` before running:

1. **Tool paths** — set `TOOLS$maple`, `TOOLS$treetime`, `TOOLS$mafft`, etc. to the correct paths on your system.
2. **Date range** — adjust `PARAMS$date_from` / `PARAMS$date_to` if needed.
3. **Clock filter threshold** — `PARAMS$clock_iqr` (default: 4 IQR) for outlier removal after round-1 TreeTime.
4. **Branch length thresholds** — `PARAMS$prune_min_bl` / `PARAMS$prune_max_bl` for MAPLE two-round pruning.

All data paths are resolved relative to the project root using `here::here()` — do not hardcode absolute paths in scripts.

---

## Running the pipeline

Run scripts in order. Each script accepts an optional `--gene` argument; omit it to run all segments.

```bash
# Stage 00 — Metadata & climate
Rscript scripts/00_metadata/00_parse_gisaid.R
Rscript scripts/00_metadata/01_geocode.R
conda activate climflu && python scripts/00_metadata/02_attach_climate.py --all

# Stage 01 — Alignment
Rscript scripts/01_alignment/01_mafft.R

# Stage 02 — Phylogenetics & dating (two-round)
Rscript scripts/02_phylogenetics/01_maple_tree.R        # round 1 tree
Rscript scripts/02_phylogenetics/02_treetime.R          # round 1 dating + clock filter
Rscript scripts/02_phylogenetics/03_maple_tree_r2.R     # round 2 tree (filtered)
Rscript scripts/02_phylogenetics/04_treetime_r2.R       # final chronogram

# Stage 03 — dN/dS estimation
Rscript scripts/03_rates_analysis/01_compute_rates.R    # overall substitution rates
Rscript scripts/03_rates_analysis/02_parsimony_dnds.R   # dN, dS, ω (whole gene)
Rscript scripts/03_rates_analysis/03_subregion_dnds.R   # HA2, epitope, ridge

# Stage 04 — Statistical analysis
Rscript scripts/03_rates_analysis/04_phylolm.R          # phylogenetic regression + FDR
Rscript scripts/03_rates_analysis/05_jackknife.R        # subsampling stability
Rscript scripts/03_rates_analysis/06_mutations.R        # recurrent mutation analysis
```

### Single gene example

```bash
Rscript scripts/00_metadata/00_parse_gisaid.R --gene H1
Rscript scripts/00_metadata/01_geocode.R
conda activate climflu && python scripts/00_metadata/02_attach_climate.py --gene H1
```

---

## Key outputs

| File | Description |
|------|-------------|
| `data/trees/{gene}_timetree.nwk` | Final time-calibrated chronogram |
| `data/rates/{gene}_rates.csv` | Per-tip overall substitution rates with climate covariates |
| `data/rates/{gene}_dnds.csv` | Per-tip dN, dS, ω and subregion rates |
| `data/rates/{gene}_phylolm_results.csv` | Phylolm coefficient tables (whole gene) |
| `data/rates/{gene}_subregion_results.csv` | Phylolm results for HA2, epitope, ridge |
| `data/rates/{gene}_jackknife.csv` | Temperature coefficients across subsamples |
| `data/rates/{gene}_mutations.csv` | Recurrent mutation–climate association tests |

---

## Acknowledgements

- Sequence data: [GISAID EpiFlu](https://gisaid.org) — please cite GISAID when using this data.
- Climate data: ERA5-Land via [Copernicus CDS](https://cds.climate.copernicus.eu). "Generated using Copernicus Climate Change Service information."
- Population density: [GHSL R2023A](https://ghsl.jrc.ec.europa.eu) (100 m, Mollweide; epochs 2010/2015/2020).
- Geocoding: [Nominatim / OpenStreetMap](https://nominatim.openstreetmap.org) via `tidygeocoder`.
- Evolutionary analysis tools: MAFFT, MAPLE, TreeTime.
- Parsimony ancestral reconstruction: phangorn R package.
- Phylogenetic regression: phylolm R package.
