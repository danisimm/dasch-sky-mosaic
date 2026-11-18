"""Builds a local, magnitude-limited, full-sky Gaia DR3 reference catalog.

Used to seed SCAMP's astrometric solve without an online catalog lookup per
plate. Queries the Gaia Archive TAP service in declination bands (to stay
under per-query row/time limits) via astroquery, then merges the bands into
a single compact FITS binary table.

Run directly (needs network access, no SExtractor/SCAMP required - can run
on either the Windows venv or the WSL venv):
    python -m dasch_sky_mosaic.gaia_catalog --mag-limit 14
"""
from __future__ import annotations

import logging
import time
from pathlib import Path

from astropy.table import Table, vstack

logging.basicConfig(level=logging.INFO, format="%(message)s")
LOG = logging.getLogger(__name__)

GAIA_TABLE = "gaiadr3.gaia_source_lite"

# gaia_source_lite has no ref_epoch/ra_error/dec_error columns. Every Gaia DR3
# source shares the same fixed reference epoch (not per-row), and the lite
# table omits position errors entirely (only proper-motion errors survive).
GAIA_DR3_REF_EPOCH = 2016.0

_COLUMNS = [
    "source_id",
    "ra",
    "dec",
    "pmra",
    "pmra_error",
    "pmdec",
    "pmdec_error",
    "phot_g_mean_mag",
]

_QUERY_RETRIES = 3
_RETRY_BACKOFF_S = 15.0


def _band_edges(dec_step_deg: float) -> list[tuple[float, float]]:
    edges = []
    dec = -90.0
    while dec < 90.0 - 1e-9:
        dec_hi = min(dec + dec_step_deg, 90.0)
        edges.append((dec, dec_hi))
        dec = dec_hi
    return edges


def _query_band(dec_min: float, dec_max: float, mag_limit: float) -> Table:
    from astroquery.gaia import Gaia

    Gaia.ROW_LIMIT = -1
    query = f"""
        SELECT {", ".join(_COLUMNS)}
        FROM {GAIA_TABLE}
        WHERE dec >= {dec_min} AND dec < {dec_max}
          AND phot_g_mean_mag < {mag_limit}
    """
    job = Gaia.launch_job_async(query)
    return job.get_results()


def download_band(
    dec_min: float,
    dec_max: float,
    mag_limit: float,
    dest_path: Path,
    overwrite: bool = False,
) -> Path:
    """Download one declination band to `dest_path`, retrying on failure."""
    if dest_path.exists() and not overwrite:
        LOG.info("Reusing cached band dec=[%.1f, %.1f): %s", dec_min, dec_max, dest_path.name)
        return dest_path

    dest_path.parent.mkdir(parents=True, exist_ok=True)
    last_exc: Exception | None = None
    for attempt in range(1, _QUERY_RETRIES + 1):
        try:
            LOG.info(
                "Querying dec=[%.1f, %.1f) G<%.1f (attempt %d/%d)",
                dec_min, dec_max, mag_limit, attempt, _QUERY_RETRIES,
            )
            table = _query_band(dec_min, dec_max, mag_limit)
            table.write(dest_path, format="fits", overwrite=True)
            LOG.info("Saved %d sources to %s", len(table), dest_path.name)
            return dest_path
        except Exception as exc:
            last_exc = exc
            LOG.warning("Band dec=[%.1f, %.1f) failed: %s", dec_min, dec_max, exc)
            if attempt < _QUERY_RETRIES:
                time.sleep(_RETRY_BACKOFF_S * attempt)

    raise RuntimeError(f"failed to download dec=[{dec_min}, {dec_max}) after {_QUERY_RETRIES} attempts") from last_exc


def merge_bands(bands_dir: Path, output_path: Path) -> Path:
    """Concatenate all per-band FITS files in `bands_dir` into one catalog."""
    band_files = sorted(bands_dir.glob("dec_*.fits"))
    if not band_files:
        raise FileNotFoundError(f"no band files found in {bands_dir}")

    LOG.info("Merging %d band files", len(band_files))
    tables = [Table.read(f) for f in band_files]
    merged = vstack(tables, metadata_conflicts="silent")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    merged.write(output_path, format="fits", overwrite=True)
    LOG.info("Wrote merged catalog: %d sources -> %s", len(merged), output_path)
    return output_path


def build_gaia_catalog(
    mag_limit: float,
    dec_step_deg: float,
    bands_dir: Path,
    output_path: Path,
    overwrite: bool = False,
) -> Path:
    """Download all declination bands (skipping cached ones) and merge them."""
    for dec_min, dec_max in _band_edges(dec_step_deg):
        band_name = f"dec_{dec_min:+06.1f}_{dec_max:+06.1f}.fits".replace("+", "p").replace("-", "m")
        download_band(dec_min, dec_max, mag_limit, bands_dir / band_name, overwrite=overwrite)

    return merge_bands(bands_dir, output_path)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Download a local, magnitude-limited Gaia DR3 catalog.")
    parser.add_argument("--mag-limit", type=float, default=14.0, help="Gaia G-band magnitude cutoff (default: 14.0)")
    parser.add_argument("--dec-step", type=float, default=10.0, help="Declination band width in degrees (default: 10.0)")
    parser.add_argument("--bands-dir", default=None, help="Directory for per-band cache files (default: data/reference_catalogs/gaia_dr3_g<mag>_bands)")
    parser.add_argument("--output", default=None, help="Output merged FITS path (default: data/reference_catalogs/gaia_dr3_g<mag>.fits)")
    parser.add_argument("--overwrite", action="store_true", help="Re-download bands even if cached")
    parser.add_argument("--merge-only", action="store_true", help="Skip downloading; just merge existing band files in --bands-dir")
    args = parser.parse_args()

    mag_tag = f"{args.mag_limit:.0f}"
    bands_dir = Path(args.bands_dir) if args.bands_dir else Path(f"data/reference_catalogs/gaia_dr3_g{mag_tag}_bands")
    output_path = Path(args.output) if args.output else Path(f"data/reference_catalogs/gaia_dr3_g{mag_tag}.fits")

    if args.merge_only:
        merge_bands(bands_dir, output_path)
    else:
        build_gaia_catalog(
            mag_limit=args.mag_limit,
            dec_step_deg=args.dec_step,
            bands_dir=bands_dir,
            output_path=output_path,
            overwrite=args.overwrite,
        )
    print(f"Done. Catalog at: {output_path}")
