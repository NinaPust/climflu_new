#!/usr/bin/env python3
"""
scripts/00_metadata/03_population_density.py

PURPOSE: Extract population density for each sample's lat/lon using the
         GHS-POP (Global Human Settlement Population Grid) rasters from JRC.
         Year is rounded to the nearest available GHS-POP epoch (every 5 years).

INPUT:
    data/metadata/{gene}_metadata.csv  (with lat/lon from 01_geocode.R)
    data/raw/pop_density/              (GHS-POP GeoTIFFs, downloaded manually)

OUTPUT:
    data/population/{gene}_popdens.csv  — metadata columns plus:
        pop_density   (people/km², from GHS-POP raster at nearest epoch)
        pop_epoch     (epoch year used, e.g. 2010, 2015, 2020)

SETUP:
    conda activate climflu
    pip install rasterio

DOWNLOAD:
    From https://ghsl.jrc.ec.europa.eu/download.php?ds=pop
    Download historical epochs (R2023A release, 100m resolution).
    Both Mollweide (54009) and WGS84 (4326) projections are supported.
    Place each downloaded folder inside data/raw/pop_density/ — the script
    searches recursively so nested folders are fine.

    For sequences collected 2009–2023, you need at minimum:
        GHS_POP_E2010_GLOBE_R2023A_54009_100_V1_0/  (or _4326_30ss_)
        GHS_POP_E2015_GLOBE_R2023A_54009_100_V1_0/
        GHS_POP_E2020_GLOBE_R2023A_54009_100_V1_0/

    NOTE: E2030 is a projected estimate, not measured historical data.
    The script will use it only if no better epoch is available.

USAGE:
    python scripts/00_metadata/03_population_density.py --gene H1
    python scripts/00_metadata/03_population_density.py --all
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
POP_RAW_DIR  = PROJECT_ROOT / "data" / "raw" / "pop_density"
POP_OUT_DIR  = PROJECT_ROOT / "data" / "population"
POP_OUT_DIR.mkdir(parents=True, exist_ok=True)

CACHE_PATH = POP_OUT_DIR / "popdens_cache.csv"

# All GHS-POP R2023A epochs (1975–2030); script uses whichever are present.
# For historical influenza data (2009–2023), download E2010, E2015, E2020.
ALL_EPOCHS = [1975, 1980, 1990, 2000, 2005, 2010, 2015, 2020, 2025, 2030]


# ── Raster helpers ────────────────────────────────────────────────────────────

def find_raster(epoch: int) -> Path | None:
    """
    Search recursively under POP_RAW_DIR for a GHS-POP raster for this epoch.
    Matches both Mollweide (54009) and WGS84 (4326) variants.
    """
    for pattern in [
        f"**/GHS_POP_E{epoch}_GLOBE_*_54009_*.tif",
        f"**/GHS_POP_E{epoch}_GLOBE_*_4326_*.tif",
        f"**/*_E{epoch}_*_54009_*.tif",
        f"**/*_E{epoch}_*_4326_*.tif",
        f"**/ghspop_{epoch}*.tif",
        f"**/pop_{epoch}*.tif",
    ]:
        hits = sorted(POP_RAW_DIR.glob(pattern))
        # Skip overview files (.ovr)
        hits = [h for h in hits if not h.suffix == ".ovr"]
        if hits:
            return hits[0]
    return None


def available_epochs() -> list[int]:
    return [e for e in ALL_EPOCHS if find_raster(e) is not None]


def round_to_epoch(year: int, epochs: list[int]) -> int:
    return min(epochs, key=lambda e: abs(e - year))


def pixel_area_km2_from_raster(src, lat: float) -> float:
    """
    Return the area in km² of one raster pixel.
    - For projected CRS (metres): area = |res_x * res_y| / 1e6
    - For geographic CRS (degrees): latitude-dependent formula
    """
    res_x, res_y = src.res  # pixel size in CRS units

    if src.crs and src.crs.is_projected:
        # Units are metres; Mollweide is equal-area so every pixel is the same
        return abs(res_x * res_y) / 1e6
    else:
        # Geographic (degrees) — area depends on latitude
        deg = abs(res_x)
        lat_rad = math.radians(lat)
        return (deg * math.pi / 180.0) ** 2 * 6371.0 ** 2 * math.cos(lat_rad)


def extract_densities_for_epoch(points_df: pd.DataFrame, epoch: int) -> pd.Series:
    """
    Extract population density (people/km²) for all points in points_df from
    the GHS-POP raster for the given epoch.  Opens the raster once and reads
    each point with a windowed single-pixel read to avoid loading the full
    global raster (~500 GB) into memory.

    points_df must have columns: lat, lon (index is preserved).
    Returns a Series of pop_density values aligned to points_df.index.
    """
    import rasterio
    from rasterio.windows import Window
    from rasterio.warp import transform as warp_transform

    tif = find_raster(epoch)
    if tif is None:
        return pd.Series(np.nan, index=points_df.index)

    results = pd.Series(np.nan, index=points_df.index)

    try:
        with rasterio.open(tif) as src:
            area_km2 = pixel_area_km2_from_raster(src, lat=0)  # equal-area: constant
            is_equal_area = src.crs and src.crs.is_projected

            for idx, row in points_df.iterrows():
                try:
                    xs, ys = warp_transform(
                        "EPSG:4326", src.crs, [row["lon"]], [row["lat"]]
                    )
                    r, c = src.index(xs[0], ys[0])
                    r = int(max(0, min(r, src.height - 1)))
                    c = int(max(0, min(c, src.width  - 1)))

                    # Read a single pixel via windowed read — avoids loading full raster
                    count = float(src.read(1, window=Window(c, r, 1, 1))[0, 0])

                    if count < 0 or np.isnan(count):
                        continue

                    if not is_equal_area:
                        area_km2 = pixel_area_km2_from_raster(src, row["lat"])

                    results[idx] = round(count / area_km2, 3) if area_km2 > 0 else np.nan

                except Exception as exc:
                    warnings.warn(
                        f"Point failed ({epoch}, {row['lat']:.3f}, "
                        f"{row['lon']:.3f}): {exc}"
                    )
    except Exception as exc:
        warnings.warn(f"Could not open raster for epoch {epoch}: {exc}")

    return results


# ── Cache management ──────────────────────────────────────────────────────────

def load_cache() -> pd.DataFrame:
    if CACHE_PATH.exists():
        cache = pd.read_csv(CACHE_PATH)
        print(f"  Loaded pop density cache: {len(cache)} entries")
        return cache
    return pd.DataFrame(columns=["lat", "lon", "pop_epoch", "pop_density"])


def save_cache(cache: pd.DataFrame):
    cache.to_csv(CACHE_PATH, index=False)


# ── Main processing ──────────────────────────────────────────────────────────

def process_genes(genes: list[str]):
    epochs = available_epochs()
    if not epochs:
        print(
            f"ERROR: No GHS-POP rasters found under {POP_RAW_DIR}.\n\n"
            "Download historical epochs from:\n"
            "  https://ghsl.jrc.ec.europa.eu/download.php?ds=pop\n\n"
            "Select Release=R2023A, Resolution=100m, Projection=Mollweide (54009)\n"
            "and download epochs E2010, E2015, E2020 (and E2000/E2005 if needed).\n"
            "Place each downloaded folder inside data/raw/pop_density/.\n\n"
            "Note: E2030 is a projected estimate, not historical data.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"Found GHS-POP rasters for epochs: {epochs}")
    historical = [e for e in epochs if e <= 2023]
    if not historical:
        print(
            "WARNING: Only future-projected epochs found (e.g. E2030). "
            "Download historical epochs (E2010, E2015, E2020) for accurate results.",
            file=sys.stderr,
        )

    gene_dfs: dict[str, pd.DataFrame] = {}
    for gene in genes:
        meta_path = METADATA_DIR / f"{gene}_metadata.csv"
        if not meta_path.exists():
            print(f"[{gene}] Metadata not found, skipping.", file=sys.stderr)
            continue
        gene_dfs[gene] = pd.read_csv(meta_path, parse_dates=["collection_date"])

    if not gene_dfs:
        print("No metadata files found. Run 00_parse_gisaid.R and "
              "01_geocode.R first.", file=sys.stderr)
        return

    all_rows = pd.concat(gene_dfs.values(), ignore_index=True)
    all_rows = all_rows.dropna(subset=["latitude", "longitude", "year"])
    all_rows["pop_epoch"] = (
        all_rows["year"].astype(int)
        .apply(lambda y: round_to_epoch(y, epochs))
    )

    unique_keys = (
        all_rows[["latitude", "longitude", "pop_epoch"]]
        .drop_duplicates()
        .rename(columns={"latitude": "lat", "longitude": "lon"})
    )
    print(f"\nTotal unique (lat, lon, epoch) combinations: {len(unique_keys)}")

    cache = load_cache()

    if not cache.empty:
        merged = unique_keys.merge(
            cache[["lat", "lon", "pop_epoch"]],
            on=["lat", "lon", "pop_epoch"],
            how="left",
            indicator=True,
        )
        to_extract = merged[merged["_merge"] == "left_only"].drop(columns="_merge")
    else:
        to_extract = unique_keys.copy()

    print(f"  Already cached: {len(unique_keys) - len(to_extract)}")
    print(f"  To extract:     {len(to_extract)}")

    if len(to_extract) > 0:
        print("\nExtracting population densities from rasters...")
        new_frames = []
        for epoch, grp in to_extract.groupby("pop_epoch"):
            print(f"  Epoch {epoch}: {len(grp)} points...")
            densities = extract_densities_for_epoch(grp.reset_index(drop=True), int(epoch))
            result = grp.reset_index(drop=True).copy()
            result["pop_density"] = densities.values
            new_frames.append(result)

        if new_frames:
            cache = pd.concat(
                [cache] + new_frames, ignore_index=True
            ).drop_duplicates(subset=["lat", "lon", "pop_epoch"], keep="last")
        save_cache(cache)
        print(f"  Cache updated: {len(cache)} total entries")

    print("\nWriting population density CSVs...")
    for gene, df in gene_dfs.items():
        out_path = POP_OUT_DIR / f"{gene}_popdens.csv"

        df = df.copy()
        df["pop_epoch"] = (
            df["year"].astype(int)
            .apply(lambda y: round_to_epoch(y, epochs))
        )
        out_df = df.merge(
            cache,
            left_on=["latitude", "longitude", "pop_epoch"],
            right_on=["lat", "lon", "pop_epoch"],
            how="left",
        ).drop(columns=["lat", "lon"], errors="ignore")

        out_df.to_csv(out_path, index=False)
        n_na = out_df["pop_density"].isna().sum()
        print(f"  [{gene}] {len(out_df)} rows, {n_na} with NA density "
              f"→ {out_path.name}")


def main():
    parser = argparse.ArgumentParser(
        description="Attach GHS-POP population density to GISAID metadata."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--gene", choices=GENES, help="Process a single gene segment"
    )
    group.add_argument(
        "--all", action="store_true", help="Process all gene segments"
    )
    args = parser.parse_args()

    process_genes(GENES if args.all else [args.gene])
    print("Done.")


if __name__ == "__main__":
    main()
