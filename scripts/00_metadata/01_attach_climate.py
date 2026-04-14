#!/usr/bin/env python3
"""
scripts/00_metadata/01_attach_climate.py

PURPOSE: Fetch ERA5-Land monthly temperature and relative humidity from
         Google Earth Engine for each sample's lat/lon/date, and attach
         to the metadata CSVs.

NOTE: This script is subject to change as metadata processing methods improve.

INPUT:
    data/metadata/{gene}_metadata.csv  (output of 00_parse_gisaid.R)

OUTPUT:
    data/climate/{gene}_climate.csv  — original metadata columns plus:
        era5_temp_mean   (°C, monthly mean 2m air temperature)
        era5_temp_min    (°C, monthly minimum)
        era5_temp_max    (°C, monthly maximum)
        era5_dewpoint    (°C, monthly mean 2m dewpoint)
        humidity_rh      (%, relative humidity from Magnus formula)

REQUIREMENTS:
    pip install earthengine-api pandas
    Run `earthengine authenticate` before first use.

USAGE:
    python scripts/00_metadata/01_attach_climate.py --gene H1
    python scripts/00_metadata/01_attach_climate.py --all
"""

import argparse
import math
import os
import sys
import warnings
from pathlib import Path

import pandas as pd

# ── Config ────────────────────────────────────────────────────────────────────
PROJECT_ROOT = Path(__file__).resolve().parents[2]
GENES = ["H1", "H3", "N1", "N2"]

METADATA_DIR = PROJECT_ROOT / "data" / "metadata"
CLIMATE_DIR  = PROJECT_ROOT / "data" / "climate"
CLIMATE_DIR.mkdir(parents=True, exist_ok=True)


def relative_humidity(temp_c: float, dewpoint_c: float) -> float:
    """
    Calculate relative humidity (%) from temperature and dewpoint
    using the Magnus formula.
    """
    a, b = 17.625, 243.04  # Magnus coefficients (°C)
    gamma_t  = (a * temp_c)    / (b + temp_c)
    gamma_dp = (a * dewpoint_c) / (b + dewpoint_c)
    rh = 100.0 * math.exp(gamma_dp - gamma_t)
    return round(rh, 2)


def fetch_era5_for_point(ee, lat: float, lon: float, year: int, month: int) -> dict:
    """
    Query ERA5-Land monthly aggregates for a single point and month.
    Returns dict with temp/dewpoint stats, or NaNs on failure.
    """
    nan_result = {
        "era5_temp_mean": float("nan"),
        "era5_temp_min":  float("nan"),
        "era5_temp_max":  float("nan"),
        "era5_dewpoint":  float("nan"),
        "humidity_rh":    float("nan"),
    }

    try:
        start = f"{year:04d}-{month:02d}-01"
        # Last day of month
        import calendar
        last_day = calendar.monthrange(year, month)[1]
        end = f"{year:04d}-{month:02d}-{last_day:02d}"

        point = ee.Geometry.Point([lon, lat])

        collection = (
            ee.ImageCollection("ECMWF/ERA5_LAND/DAILY_AGGR")
            .filterDate(start, end)
            .filterBounds(point)
            .select(["temperature_2m", "dewpoint_temperature_2m"])
        )

        # Monthly stats
        temp_mean = collection.select("temperature_2m").mean()
        temp_min  = collection.select("temperature_2m").min()
        temp_max  = collection.select("temperature_2m").max()
        dew_mean  = collection.select("dewpoint_temperature_2m").mean()

        def sample_val(img):
            val = img.reduceRegion(
                reducer=ee.Reducer.first(),
                geometry=point,
                scale=11132,  # ERA5-Land native resolution ~0.1°
            )
            return val.getInfo()

        t_mean = sample_val(temp_mean).get("temperature_2m")
        t_min  = sample_val(temp_min).get("temperature_2m")
        t_max  = sample_val(temp_max).get("temperature_2m")
        dp     = sample_val(dew_mean).get("dewpoint_temperature_2m")

        if any(v is None for v in [t_mean, t_min, t_max, dp]):
            return nan_result

        # ERA5 values are in Kelvin — convert to Celsius
        t_mean_c = t_mean - 273.15
        t_min_c  = t_min  - 273.15
        t_max_c  = t_max  - 273.15
        dp_c     = dp     - 273.15
        rh       = relative_humidity(t_mean_c, dp_c)

        return {
            "era5_temp_mean": round(t_mean_c, 3),
            "era5_temp_min":  round(t_min_c, 3),
            "era5_temp_max":  round(t_max_c, 3),
            "era5_dewpoint":  round(dp_c, 3),
            "humidity_rh":    rh,
        }

    except Exception as exc:
        warnings.warn(f"GEE query failed for ({lat},{lon}) {year}-{month}: {exc}")
        return nan_result


def process_gene(ee, gene: str):
    meta_path = METADATA_DIR / f"{gene}_metadata.csv"
    if not meta_path.exists():
        print(f"[{gene}] Metadata not found, skipping: {meta_path}", file=sys.stderr)
        return

    out_path = CLIMATE_DIR / f"{gene}_climate.csv"
    print(f"[{gene}] Processing {meta_path}...")

    df = pd.read_csv(meta_path, parse_dates=["collection_date"])
    df["year"]  = df["collection_date"].dt.year
    df["month"] = df["collection_date"].dt.month

    climate_rows = []
    for _, row in df.iterrows():
        result = fetch_era5_for_point(
            ee,
            lat=row["latitude"],
            lon=row["longitude"],
            year=int(row["year"]),
            month=int(row["month"]),
        )
        climate_rows.append(result)

    climate_df = pd.DataFrame(climate_rows)
    out_df = pd.concat([df.reset_index(drop=True), climate_df], axis=1)
    out_df.to_csv(out_path, index=False)
    print(f"[{gene}] Written: {out_path}  ({len(out_df)} rows)")


def main():
    parser = argparse.ArgumentParser(description="Attach ERA5 climate data to GISAID metadata.")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gene", choices=GENES, help="Process a single gene segment")
    group.add_argument("--all", action="store_true", help="Process all gene segments")
    args = parser.parse_args()

    # Authenticate and initialise GEE
    try:
        import ee
        ee.Initialize()
    except Exception as exc:
        sys.exit(
            f"Google Earth Engine initialisation failed: {exc}\n"
            "Run `earthengine authenticate` and retry."
        )

    genes_to_process = GENES if args.all else [args.gene]
    for gene in genes_to_process:
        process_gene(ee, gene)

    print("Done.")


if __name__ == "__main__":
    main()
