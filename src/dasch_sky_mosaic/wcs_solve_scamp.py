"""Solve WCS for plate JPGs using SExtractor + SCAMP against a local Gaia catalog.

Alternative to wcs_solve.py's astroalign backend: needs only the plate JPG (no
FITS mosaic download), at the cost of a rougher initial guess that sometimes
fails to solve (see scamp_prep.py for where that guess comes from and its
known limits).

Run directly:
    python -m dasch_sky_mosaic.wcs_solve_scamp plate_id1 plate_id2 \\
        --exposures dev/exposures_single.csv \\
        --ref-catalog data/reference_catalogs/gaia_dr3_g14_ldac.fits

Only runs where SExtractor/SCAMP are installed (WSL) - see README.
"""
from __future__ import annotations

import csv
import logging
import math
import os
import subprocess
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.io.votable import parse as parse_votable
from PIL import Image

from dasch_sky_mosaic.call_s3 import download_photo_from_s3
from dasch_sky_mosaic.scamp_prep import initial_guess_wcs, jpg_to_fits

logging.basicConfig(level=logging.INFO, format="%(message)s")
LOG = logging.getLogger(__name__)

# Relative, not absolute: SExtractor/SCAMP mishandle spaces in absolute paths
# (this repo lives under a directory with a space in it), so config file
# arguments must be given relative to the current working directory instead.
# Assumes the script is run from the repo root, like the other CLI scripts.
CONFIG_DIR = Path(os.path.relpath(Path(__file__).parent / "scamp_config", Path.cwd()))
ROTATION_HYPOTHESES_DEG = (0.0, 90.0, 180.0, 270.0)

# Must stay >= scamp.conf's POSITION_MAXERR (currently 360 arcmin = 6 deg) -
# this bounds how far the reference catalog crop reaches, and cropping
# tighter than SCAMP's own search radius would silently exclude valid stars.
POSITION_MAXERR_DEG = 6.0

# Reduced chi2 (chi2 / n_matched) gets statistically meaningless at very low
# n - a "solution" with a handful of coincidental matches can show a low raw
# chi2 despite being a false positive (see the b02312 comparison: the 3 wrong
# rotations each matched ~8 detections vs. the true one's 489). Below this
# floor a hypothesis is rejected outright regardless of its chi2.
MIN_MATCHED_DETECTIONS = 20


def _rotate_cd(header: dict[str, Any], rot_deg: float) -> dict[str, Any]:
    """Return a copy of `header` with its CD matrix rotated by `rot_deg`.

    CRVAL/CRPIX are untouched - only the assumed orientation changes.
    """
    cd = np.array([[header["CD1_1"], header["CD1_2"]], [header["CD2_1"], header["CD2_2"]]])
    rot = math.radians(rot_deg)
    r = np.array([[math.cos(rot), -math.sin(rot)], [math.sin(rot), math.cos(rot)]])
    cd_rot = cd @ r

    rotated = dict(header)
    rotated["CD1_1"], rotated["CD1_2"] = float(cd_rot[0, 0]), float(cd_rot[0, 1])
    rotated["CD2_1"], rotated["CD2_2"] = float(cd_rot[1, 0]), float(cd_rot[1, 1])
    return rotated


def crop_reference_catalog(
    ref_catalog_path: Path,
    center_ra_deg: float,
    center_dec_deg: float,
    radius_deg: float,
    out_path: Path,
) -> Path:
    """Write a local subset of the reference catalog within `radius_deg` of the guess center.

    SCAMP has to search the whole reference catalog for every detection - our
    full-sky catalog has ~17 million stars, most of them nowhere near any
    given plate. Cropping to just the region a plate could plausibly cover
    (position guess error + field size) is the dominant lever on solve time.
    """
    with fits.open(ref_catalog_path) as hdul:
        data = hdul["LDAC_OBJECTS"].data

        ra1, dec1 = math.radians(center_ra_deg), math.radians(center_dec_deg)
        ra2, dec2 = np.radians(data["X_WORLD"]), np.radians(data["Y_WORLD"])
        cos_sep = np.sin(dec1) * np.sin(dec2) + np.cos(dec1) * np.cos(dec2) * np.cos(ra1 - ra2)
        sep_deg = np.degrees(np.arccos(np.clip(cos_sep, -1.0, 1.0)))

        cropped = data[sep_deg <= radius_deg]
        columns = [
            fits.Column(name=col.name, format=col.format, unit=col.unit, array=cropped[col.name])
            for col in data.columns
        ]
        fits.HDUList([fits.PrimaryHDU(), fits.BinTableHDU.from_columns(columns, name="LDAC_OBJECTS")]).writeto(
            out_path, overwrite=True
        )

    return out_path


def run_sextractor(fits_path: Path, out_ldac_path: Path) -> None:
    out_ldac_path.parent.mkdir(parents=True, exist_ok=True)
    result = subprocess.run(
        [
            "source-extractor", str(fits_path),
            "-c", str(CONFIG_DIR / "sextractor.sex"),
            "-PARAMETERS_NAME", str(CONFIG_DIR / "sextractor.param"),
            "-FILTER_NAME", str(CONFIG_DIR / "sextractor.conv"),
            "-CATALOG_NAME", str(out_ldac_path),
        ],
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        raise RuntimeError(f"source-extractor failed for {fits_path.name}:\n{result.stderr}")


def _parse_scamp_xml(xml_path: Path) -> tuple[int, float] | None:
    """Return (n_matched_detections, chi2) against the reference catalog from SCAMP's XML report."""
    if not xml_path.exists():
        return None
    votable = parse_votable(xml_path)
    for table in votable.iter_tables():
        if table.name == "FGroups":
            row = table.array[0]
            return int(row["AstromNDets_Reference"]), float(row["AstromChi2_Reference"])
    return None


def run_scamp(ldac_path: Path, ref_catalog_path: Path) -> tuple[Path, int, float] | None:
    """Run SCAMP on `ldac_path`. Returns (head_path, n_matched_detections, chi2), or None on failure."""
    head_path = ldac_path.with_suffix(".head")
    xml_path = ldac_path.with_suffix(".xml")
    result = subprocess.run(
        [
            "scamp", str(ldac_path),
            "-c", str(CONFIG_DIR / "scamp.conf"),
            "-ASTREFCAT_NAME", str(ref_catalog_path),
            "-XML_NAME", str(xml_path),
        ],
        capture_output=True,
        text=True,
    )
    # SCAMP's own match-count/chi2 reporting goes to stderr regardless of
    # outcome - always logging it is how we judge whether a "success" is a
    # trustworthy fit rather than a false positive.
    LOG.info("SCAMP output for %s:\n%s", ldac_path.name, result.stderr.strip()[-2000:])

    if result.returncode != 0:
        LOG.info("SCAMP exit code for %s: %d", ldac_path.name, result.returncode)
        return None
    if not head_path.exists():
        LOG.info("SCAMP produced no .head for %s", ldac_path.name)
        return None

    stats = _parse_scamp_xml(xml_path)
    if stats is None:
        LOG.info("SCAMP produced no parseable XML stats for %s", ldac_path.name)
        return None
    return head_path, stats[0], stats[1]


def _header_to_jsonable(header: fits.Header) -> dict[str, Any]:
    """Convert a FITS header to a plain JSON-safe dict.

    dict(header) breaks on repeated COMMENT/HISTORY cards (astropy collapses
    them into a non-serializable _HeaderCommentaryCards object) - iterating
    header.items() instead yields each one individually, so they're gathered
    into a list per key rather than dropped. SCAMP often reports fit-quality
    info (matches, chi2) as COMMENT lines, useful for comparing hypotheses.
    """
    result: dict[str, Any] = {}
    for key, value in header.items():
        if key == "":
            continue
        if key in ("COMMENT", "HISTORY"):
            result.setdefault(key, []).append(str(value))
        else:
            result[key] = value
    return result


def solve_plate_wcs_scamp(
    plate_id: str,
    jpg_path: Path,
    exposures_by_plate: dict[str, dict[str, str]],
    ref_catalog_path: Path,
    work_dir: Path,
) -> dict[str, Any]:
    plate_dir = work_dir / plate_id
    plate_dir.mkdir(parents=True, exist_ok=True)

    with Image.open(jpg_path) as im:
        jpg_w, jpg_h = im.size

    base_header = initial_guess_wcs(plate_id, jpg_w, jpg_h, exposures_by_plate)

    # Crop once (same center for all 4 rotations) - the dominant cost in a
    # SCAMP run is searching the reference catalog, so this alone is roughly
    # a 200x reduction in candidate stars vs. the full 17-million-star catalog.
    scale_deg_px = math.sqrt(abs(np.linalg.det(np.array([
        [base_header["CD1_1"], base_header["CD1_2"]],
        [base_header["CD2_1"], base_header["CD2_2"]],
    ]))))
    field_radius_deg = 0.5 * scale_deg_px * math.hypot(jpg_w, jpg_h)
    crop_radius_deg = POSITION_MAXERR_DEG + field_radius_deg + 0.5

    cropped_ref_catalog = plate_dir / "ref_catalog_crop.fits"
    crop_reference_catalog(
        ref_catalog_path, base_header["CRVAL1"], base_header["CRVAL2"], crop_radius_deg, cropped_ref_catalog
    )

    all_hypotheses = []
    for rot_deg in ROTATION_HYPOTHESES_DEG:
        hypothesis_header = _rotate_cd(base_header, rot_deg)
        hypothesis_fits = plate_dir / f"plate_rot{int(rot_deg)}.fits"
        hypothesis_ldac = plate_dir / f"plate_rot{int(rot_deg)}.ldac"

        jpg_to_fits(jpg_path, hypothesis_fits, header=hypothesis_header)
        run_sextractor(hypothesis_fits, hypothesis_ldac)

        scamp_result = run_scamp(hypothesis_ldac, cropped_ref_catalog)
        if scamp_result is None:
            LOG.info("No solution for %s at rotation hypothesis %.0f deg", plate_id, rot_deg)
            all_hypotheses.append({"rotation_hypothesis_deg": rot_deg, "solved": False})
            continue

        head_path, n_matched, chi2 = scamp_result
        reduced_chi2 = chi2 / n_matched if n_matched else float("inf")
        solved_header = fits.Header.fromtextfile(head_path)
        wcs_header = _header_to_jsonable(solved_header)
        LOG.info(
            "Solved %s at rotation hypothesis %.0f deg (n_matched=%d, chi2=%.1f, reduced_chi2=%.3f)",
            plate_id, rot_deg, n_matched, chi2, reduced_chi2,
        )
        all_hypotheses.append({
            "rotation_hypothesis_deg": rot_deg,
            "solved": True,
            "n_matched": n_matched,
            "chi2": chi2,
            "reduced_chi2": reduced_chi2,
            "wcs_header": wcs_header,
        })

    candidates = [h for h in all_hypotheses if h["solved"] and h["n_matched"] >= MIN_MATCHED_DETECTIONS]
    if not candidates:
        raise RuntimeError(
            f"SCAMP found no solution for {plate_id} with >= {MIN_MATCHED_DETECTIONS} matched detections "
            "across all rotation hypotheses"
        )

    if len(candidates) > 1:
        LOG.warning(
            "%s solved at %d different rotation hypotheses (%s) - picking lowest reduced chi2",
            plate_id, len(candidates), [h["rotation_hypothesis_deg"] for h in candidates],
        )

    best = min(candidates, key=lambda h: h["reduced_chi2"])
    return {
        "plate_id": plate_id,
        "image_width_px": jpg_w,
        "image_height_px": jpg_h,
        "wcs_header": best["wcs_header"],
        "alignment_meta": {
            "backend": "scamp",
            "rotation_hypothesis_deg": best["rotation_hypothesis_deg"],
            "n_matched": best["n_matched"],
            "reduced_chi2": best["reduced_chi2"],
            "n_hypotheses_solved": len(candidates),
        },
        "all_hypotheses": all_hypotheses,
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description="Solve WCS for DASCH plates via SExtractor + SCAMP.")
    parser.add_argument("plate_ids", nargs="+", help="Plate IDs to process (e.g. b02312 a11740)")
    parser.add_argument("--exposures", default="dev/exposures_single.csv", help="Path to the cleaned single-exposure catalog")
    parser.add_argument(
        "--ref-catalog",
        default="data/reference_catalogs/gaia_dr3_g14_ldac.fits",
        help="Path to the local Gaia LDAC reference catalog",
    )
    parser.add_argument("--photo-dir", default="data/cache/photos", help="Directory to cache downloaded JPGs")
    parser.add_argument("--work-dir", default="data/cache/scamp_solve", help="Working directory for intermediate files")
    parser.add_argument("--output-dir", default="data/output/wcs_scamp", help="Directory for output JSON files")
    args = parser.parse_args()

    with open(args.exposures, newline="", encoding="utf-8") as f:
        exposures_by_plate = {row["plateid"]: row for row in csv.DictReader(f)}

    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    solved = 0
    for plate_id in args.plate_ids:
        plate_id = plate_id.strip().lower()
        LOG.info("Processing %s", plate_id)
        try:
            jpg_path = download_photo_from_s3(plate_id, Path(args.photo_dir))
            result = solve_plate_wcs_scamp(
                plate_id, jpg_path, exposures_by_plate, Path(args.ref_catalog), Path(args.work_dir)
            )
            (output_dir / f"{plate_id}_wcs.json").write_text(json.dumps(result, indent=2))
            solved += 1
        except Exception as exc:
            LOG.error("Failed for %s: %s", plate_id, exc)

    print(f"Done. {solved}/{len(args.plate_ids)} plates solved.")
