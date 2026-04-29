#!/usr/bin/env python3
"""
scripts/00_metadata/02_attach_climate.py

PURPOSE: Fetch ERA5-Land monthly temperature and relative humidity from
         the Copernicus Climate Data Store (CDS) for each sample's
         lat/lon/date, and attach to the metadata CSVs.

         Downloads one NetCDF per unique year-month (bounding box of all
         points in that month), extracts point values, and caches results.
         Re-runs only fetch new months.

INPUT:
    data/metadata/{gene}_metadata.csv  (with lat/lon from 01_geocode.R)

OUTPUT:
    data/climate/{gene}_climate.csv  — original metadata columns plus:
        era5_temp_mean   (°C, monthly mean 2m air temperature)
        era5_dewpoint    (°C, monthly mean 2m dewpoint)
        humidity_rh      (%, relative humidity from Magnus formula)

SETUP:
    1. Register at https://cds.climate.copernicus.eu
    2. Accept licence at the ERA5-Land monthly means dataset page
    3. Create ~/.cdsapirc with your key (see CDS docs)

    conda activate climflu
    pip install cdsapi pandas netcdf4 xarray

USAGE:
    python scripts/00_metadata/02_attach_climate.py --gene H1
    python scripts/00_metadata/02_attach_climate.py --all
"""

import argparse
import math
import sys
import warnings
from pathlib import Path

import pandas as pd
import numpy as np

# ── Config ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
GENES = ["H1", "H3", "N1", "N2", "PA_H1N1", "PA_H3N2"]

METADATA_DIR = PROJECT_ROOT / "data" / "metadata"
CLIMATE_DIR = PROJECT_ROOT / "data" / "climate"
NC_DIR = CLIMATE_DIR / "nc_downloads"
CLIMATE_DIR.mkdir(parents=True, exist_ok=True)
NC_DIR.mkdir(parents=True, exist_ok=True)

CACHE_PATH = CLIMATE_DIR / "era5_cache.csv"

CLIMATE_COLS = [
    "era5_temp_mean",
    "era5_dewpoint",
    "humidity_rh",
]

# Columns stored for the previous month (same structure, different prefix)
PREV_CLIMATE_COLS = [
    "era5_temp_prev",
    "era5_dewpoint_prev",
    "humidity_rh_prev",
]


def relative_humidity(temp_c: float, dewpoint_c: float) -> float:
    """
    Calculate relative humidity (%) from temperature and dewpoint
    using the Magnus formula.
    """
    a, b = 17.625, 243.04
    gamma_t = (a * temp_c) / (b + temp_c)
    gamma_dp = (a * dewpoint_c) / (b + dewpoint_c)
    rh = 100.0 * math.exp(gamma_dp - gamma_t)
    return round(rh, 2)


# Vectorised version for arrays
def relative_humidity_vec(temp_c, dewpoint_c):
    a, b = 17.625, 243.04
    gamma_t = (a * temp_c) / (b + temp_c)
    gamma_dp = (a * dewpoint_c) / (b + dewpoint_c)
    return np.round(100.0 * np.exp(gamma_dp - gamma_t), 2)


def download_era5_month(client, year: int, month: int,
                        lat_min: float, lat_max: float,
                        lon_min: float, lon_max: float) -> Path:
    """Download ERA5-Land monthly means for one month, return path to NetCDF."""
    nc_file = NC_DIR / f"era5land_{year:04d}_{month:02d}.nc"

    if nc_file.exists():
        return nc_file

    # Pad bounding box by 1°
    area = [
        min(lat_max + 1, 90),
        max(lon_min - 1, -180),
        max(lat_min - 1, -90),
        min(lon_max + 1, 180),
    ]

    print(f"    Downloading ERA5-Land for {year}-{month:02d} "
          f"(area: {area})...")

    zip_file = nc_file.with_suffix(".zip")

    client.retrieve(
        "reanalysis-era5-land-monthly-means",
        {
            "product_type": ["monthly_averaged_reanalysis"],
            "variable": ["2m_temperature", "2m_dewpoint_temperature"],
            "year": [str(year)],
            "month": [f"{month:02d}"],
            "time": ["00:00"],
            "area": area,
            "data_format": "netcdf",
        },
        str(zip_file),
    )

    # CDS delivers NetCDF inside a zip archive — extract it
    import zipfile
    with zipfile.ZipFile(zip_file, "r") as zf:
        nc_names = [n for n in zf.namelist() if n.endswith(".nc")]
        if not nc_names:
            raise RuntimeError(f"No .nc file found in {zip_file}")
        zf.extract(nc_names[0], NC_DIR)
        extracted = NC_DIR / nc_names[0]
        if extracted != nc_file:
            extracted.rename(nc_file)
    zip_file.unlink()

    return nc_file


def extract_points_from_nc(nc_file: Path, points_df: pd.DataFrame) -> pd.DataFrame:
    """Extract temperature and dewpoint at each (lat, lon) from a NetCDF file."""
    import xarray as xr

    ds = xr.open_dataset(nc_file)

    # ERA5-Land NetCDF variable names
    # Try common names — CDS uses t2m/d2m or full names
    temp_var = None
    dew_var = None
    for name in ds.data_vars:
        lower = name.lower()
        if "dewpoint" in lower or name == "d2m":
            dew_var = name
        elif "temperature" in lower or name == "t2m":
            temp_var = name

    if temp_var is None or dew_var is None:
        warnings.warn(f"Cannot identify variables in {nc_file}. "
                      f"Found: {list(ds.data_vars)}")
        points_df["era5_temp_mean"] = np.nan
        points_df["era5_dewpoint"] = np.nan
        points_df["humidity_rh"] = np.nan
        return points_df

    temp_vals = []
    dew_vals = []
    for _, row in points_df.iterrows():
        # Select nearest grid point
        try:
            t = float(ds[temp_var].sel(
                latitude=row["lat"], longitude=row["lon"],
                method="nearest"
            ).values.flat[0])
            d = float(ds[dew_var].sel(
                latitude=row["lat"], longitude=row["lon"],
                method="nearest"
            ).values.flat[0])
        except Exception:
            t, d = np.nan, np.nan
        temp_vals.append(t)
        dew_vals.append(d)

    ds.close()

    temp_c = np.array(temp_vals) - 273.15  # Kelvin → Celsius
    dew_c = np.array(dew_vals) - 273.15

    points_df = points_df.copy()
    points_df["era5_temp_mean"] = np.round(temp_c, 3)
    points_df["era5_dewpoint"] = np.round(dew_c, 3)
    points_df["humidity_rh"] = relative_humidity_vec(temp_c, dew_c)

    return points_df


# ── Cache management ──────────────────────────────────────────────────────────

def prev_year_month(year: int, month: int) -> tuple[int, int]:
    """Return (year, month) for the calendar month before the given one."""
    if month == 1:
        return year - 1, 12
    return year, month - 1


def load_cache() -> pd.DataFrame:
    if CACHE_PATH.exists():
        cache = pd.read_csv(CACHE_PATH)
        print(f"  Loaded ERA5 cache: {len(cache)} entries")
        return cache
    return pd.DataFrame(columns=["lat", "lon", "year", "month"] + CLIMATE_COLS)


def save_cache(cache: pd.DataFrame):
    cache.to_csv(CACHE_PATH, index=False)


# ── Main processing ──────────────────────────────────────────────────────────

def process_genes(genes: list[str]):
    import cdsapi
    client = cdsapi.Client()

    # Load all requested metadata
    gene_dfs = {}
    for gene in genes:
        meta_path = METADATA_DIR / f"{gene}_metadata.csv"
        if not meta_path.exists():
            print(f"[{gene}] Metadata not found, skipping: {meta_path}",
                  file=sys.stderr)
            continue
        df = pd.read_csv(meta_path, parse_dates=["collection_date"])
        gene_dfs[gene] = df

    if not gene_dfs:
        print("No metadata files found. Run 00_parse_gisaid.R and "
              "01_geocode.R first.", file=sys.stderr)
        return

    # Collect all (lat, lon, year, month) combos across genes, plus previous months
    all_rows = pd.concat(gene_dfs.values(), ignore_index=True)
    all_rows = all_rows.dropna(subset=["latitude", "longitude", "year", "month"])

    unique_keys = (
        all_rows[["latitude", "longitude", "year", "month"]]
        .drop_duplicates()
        .rename(columns={"latitude": "lat", "longitude": "lon"})
    )
    unique_keys["year"] = unique_keys["year"].astype(int)
    unique_keys["month"] = unique_keys["month"].astype(int)

    # Expand to also include each location's previous calendar month
    prev_keys = unique_keys.copy()
    prev_ym = [prev_year_month(r.year, r.month)
               for r in prev_keys.itertuples(index=False)]
    prev_keys["year"]  = [ym[0] for ym in prev_ym]
    prev_keys["month"] = [ym[1] for ym in prev_ym]

    all_unique_keys = (
        pd.concat([unique_keys, prev_keys], ignore_index=True)
        .drop_duplicates()
    )

    print(f"\nTotal unique (lat, lon, year, month) combinations "
          f"(current + previous months): {len(all_unique_keys)}")

    # Load cache and find uncached keys
    cache = load_cache()
    if not cache.empty:
        merged = all_unique_keys.merge(
            cache[["lat", "lon", "year", "month"]],
            on=["lat", "lon", "year", "month"],
            how="left",
            indicator=True,
        )
        to_fetch = merged[merged["_merge"] == "left_only"].drop(
            columns="_merge"
        )
    else:
        to_fetch = all_unique_keys

    print(f"  Already cached: {len(unique_keys) - len(to_fetch)}")
    print(f"  To fetch: {len(to_fetch)}")

    # Group by year-month and process
    if len(to_fetch) > 0:
        grouped = to_fetch.groupby(["year", "month"])
        n_months = len(grouped)
        print(f"\nDownloading ERA5-Land for {n_months} year-months...")

        all_new = []
        for i, ((yr, mn), grp) in enumerate(grouped, 1):
            print(f"  [{i}/{n_months}] {yr}-{mn:02d} "
                  f"({len(grp)} points)")

            try:
                nc_file = download_era5_month(
                    client, yr, mn,
                    grp["lat"].min(), grp["lat"].max(),
                    grp["lon"].min(), grp["lon"].max(),
                )
                result = extract_points_from_nc(nc_file, grp)
            except Exception as exc:
                warnings.warn(f"Failed for {yr}-{mn}: {exc}")
                result = grp.copy()
                for col in CLIMATE_COLS:
                    result[col] = np.nan

            all_new.append(result)

            # Checkpoint every 10 months
            if i % 10 == 0:
                cache = pd.concat([cache] + all_new, ignore_index=True)
                cache = cache.drop_duplicates(
                    subset=["lat", "lon", "year", "month"], keep="last"
                )
                save_cache(cache)
                all_new = []
                print(f"    Cache checkpoint: {len(cache)} entries")

        if all_new:
            cache = pd.concat([cache] + all_new, ignore_index=True)
            cache = cache.drop_duplicates(
                subset=["lat", "lon", "year", "month"], keep="last"
            )
        save_cache(cache)
        print(f"  Cache updated: {len(cache)} total entries")

    # Join climate data back to each gene's metadata
    print("\nWriting climate CSVs...")
    for gene, df in gene_dfs.items():
        out_path = CLIMATE_DIR / f"{gene}_climate.csv"

        # Current-month climate
        out_df = df.merge(
            cache,
            left_on=["latitude", "longitude", "year", "month"],
            right_on=["lat", "lon", "year", "month"],
            how="left",
        ).drop(columns=["lat", "lon"], errors="ignore")

        # Compute previous (year, month) for each row
        prev_ym = [prev_year_month(int(r.year), int(r.month))
                   for r in out_df[["year", "month"]].itertuples(index=False)]
        out_df["_prev_year"]  = [ym[0] for ym in prev_ym]
        out_df["_prev_month"] = [ym[1] for ym in prev_ym]

        # Previous-month climate
        prev_cache = cache.rename(columns={
            "year":          "_prev_year",
            "month":         "_prev_month",
            "era5_temp_mean": "era5_temp_prev",
            "era5_dewpoint":  "era5_dewpoint_prev",
            "humidity_rh":    "humidity_rh_prev",
        })

        out_df = out_df.merge(
            prev_cache[["lat", "lon", "_prev_year", "_prev_month",
                         "era5_temp_prev", "era5_dewpoint_prev", "humidity_rh_prev"]],
            left_on=["latitude", "longitude", "_prev_year", "_prev_month"],
            right_on=["lat", "lon", "_prev_year", "_prev_month"],
            how="left",
        ).drop(columns=["lat", "lon", "_prev_year", "_prev_month"], errors="ignore")

        # Differences: current minus previous month
        out_df["delta_temp"]     = np.round(
            out_df["era5_temp_mean"] - out_df["era5_temp_prev"], 3)
        out_df["delta_humidity"] = np.round(
            out_df["humidity_rh"] - out_df["humidity_rh_prev"], 3)

        n_before = len(out_df)
        out_df = out_df.dropna(subset=["era5_temp_mean"])
        out_df.to_csv(out_path, index=False)
        n_no_prev = out_df["era5_temp_prev"].isna().sum()
        print(f"  [{gene}] {len(out_df)}/{n_before} rows kept "
              f"(dropped {n_before - len(out_df)} with NA climate); "
              f"{n_no_prev} rows with NA previous-month climate "
              f"→ {out_path.name}")


def main():
    parser = argparse.ArgumentParser(
        description="Attach ERA5 climate data to GISAID metadata."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--gene", choices=GENES, help="Process a single gene segment"
    )
    group.add_argument(
        "--all", action="store_true", help="Process all gene segments"
    )
    args = parser.parse_args()

    genes_to_process = GENES if args.all else [args.gene]
    process_genes(genes_to_process)

    print("Done.")


if __name__ == "__main__":
    main()
