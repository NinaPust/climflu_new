# climflu_new

Studying the effects of environmental variables (temperature, relative humidity) on the molecular evolution of Influenza A surface proteins — a redo of [Influenza_env](https://github.com/NinaPust/Influenza_env) with GISAID data and improved methods.

**Gene segments:** H1, H3, N1, N2  
**Environmental data:** ERA5-Land (Google Earth Engine)  
**Statistical approach:** phylogenetic comparative regression (phylolm)

---

## Pipeline overview

```
Stage 00 — Metadata
  GISAID export (manual) → parse metadata → attach ERA5 climate

Stage 01 — Alignment
  MAFFT multiple alignment → MACSE pairwise codon alignment

Stage 02 — Phylogenetics
  MAPLE tree → codeml dN/dS → assemble distance matrices → TreeTime dating

Stage 03 — Rates & analysis
  Extract substitution rates → phylolm (rates & dN/dS ~ temp + humidity)
```

---

## Data

**GISAID sequences cannot be downloaded programmatically.** You must export them manually from [gisaid.org](https://gisaid.org).

### Required exports (one per gene segment: H1, H3, N1, N2)

From the EpiFlu database, download:

| File | Naming convention |
|------|-------------------|
| Sequence FASTA | `data/raw/{gene}_sequences.fasta` |
| Metadata TSV | `data/raw/{gene}_metadata.tsv` |

**Metadata fields to include when downloading** (select in the GISAID export dialog):  
Accession ID, Virus name, Collection date, Location, Latitude, Longitude, Host

These files are gitignored and must be obtained independently.

---

## Requirements

### R packages

```r
install.packages(c("here", "dplyr", "readr", "stringr", "lubridate",
                   "ape", "phangorn", "phylolm", "MuMIn"))
```

Or restore the full environment:

```r
# install.packages("renv")
renv::restore()
```

### Python packages

```bash
pip install earthengine-api pandas
earthengine authenticate   # first time only
```

### External CLI tools

| Tool | Purpose | Version tested |
|------|---------|----------------|
| [MAFFT](https://mafft.cbrc.jp/alignment/software/) | Multiple sequence alignment | ≥ 7.5 |
| [MACSE](https://bioweb.supagro.inra.fr/macse/) | Codon-aware alignment | v2.07 |
| [MAPLE](https://github.com/NicolaDM/MAPLE) | Phylogenetic inference | latest |
| [PAML](http://abacus.gene.ucl.ac.uk/software/paml.html) (codeml) | dN/dS estimation | 4.8+ |
| [TreeTime](https://github.com/neherlab/treetime) | Molecular clock dating | ≥ 0.9 |

---

## Configuration

Edit `config.R` before running:

1. **Tool paths** — set `TOOLS$macse`, `TOOLS$maple`, etc. to the correct paths on your system.
2. **Date range** — adjust `PARAMS$date_from` / `PARAMS$date_to` if needed.
3. **Branch length thresholds** — `PARAMS$prune_min_bl` / `PARAMS$prune_max_bl` for MAPLE pruning.

All data paths are resolved relative to the project root using `here::here()` — do not hardcode absolute paths in scripts.

---

## Running the pipeline

Run scripts in order. Each script accepts an optional `--gene` argument; omit it to run all four segments.

```bash
# Stage 00 — Metadata
Rscript scripts/00_metadata/00_parse_gisaid.R
python  scripts/00_metadata/01_attach_climate.py --all

# Stage 01 — Alignment
Rscript scripts/01_alignment/01_mafft.R
Rscript scripts/02_phylogenetics/01_maple_tree.R   # needed before macse
Rscript scripts/01_alignment/02_macse.R

# Stage 02 — Phylogenetics
Rscript scripts/02_phylogenetics/02_codeml.R
Rscript scripts/02_phylogenetics/03_assemble_matrices.R
Rscript scripts/02_phylogenetics/04_treetime.R

# Stage 03 — Rates & analysis
Rscript scripts/03_rates_analysis/01_compute_rates.R
Rscript scripts/03_rates_analysis/02_phylolm.R
```

### Single gene example

```bash
Rscript scripts/00_metadata/00_parse_gisaid.R --gene H1
python  scripts/00_metadata/01_attach_climate.py --gene H1
```

---

## Key outputs

| File | Description |
|------|-------------|
| `data/trees/{gene}_timetree.nwk` | Time-calibrated molecular clock tree |
| `data/trees/{gene}_dN_tree.nwk` | dN distance tree |
| `data/trees/{gene}_dS_tree.nwk` | dS distance tree |
| `data/rates/{gene}_rates.csv` | Per-tip rates, dN/dS, and environmental variables |
| `data/rates/{gene}_phylolm_results.csv` | Phylolm coefficient tables |

---

## Acknowledgements

- Sequence data: [GISAID EpiFlu](https://gisaid.org) — please cite GISAID when using this data.
- Climate data: ERA5-Land via [Google Earth Engine](https://earthengine.google.com).
- Evolutionary analysis tools: MAFFT, MACSE, MAPLE, PAML, TreeTime.
