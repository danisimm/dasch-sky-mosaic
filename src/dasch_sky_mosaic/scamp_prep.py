"""Rough WCS guess and JPG->FITS conversion, to prepare a plate for SExtractor + SCAMP.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits
from PIL import Image

# Empirically derived.
_SCALE_CORRECTION = 1.19


def initial_guess_wcs(
    plate_id: str,
    jpg_width_px: int,
    jpg_height_px: int,
    exposures_by_plate: dict[str, dict[str, str]],
) -> dict[str, Any]:
    """Return a rough initial-guess WCS header dict for `plate_id`.
    """
    row = exposures_by_plate[plate_id]

    delta_ra_x = float(row["delta_ra_x"])
    delta_ra_y = float(row["delta_ra_y"])
    delta_dec_x = float(row["delta_dec_x"])
    delta_dec_y = float(row["delta_dec_y"])
    naxis1 = float(row["naxis1"])
    naxis2 = float(row["naxis2"])

    cd_raw = np.array([[delta_ra_x, delta_ra_y], [delta_dec_x, delta_dec_y]])
    raw_det = np.linalg.det(cd_raw)

    scale_deg_px = math.sqrt(abs(raw_det) * naxis1 * naxis2 / (jpg_width_px * jpg_height_px)) * _SCALE_CORRECTION

    # Every sampled plate needed a mirror flip relative to the raw grid.
    theta = math.atan2(delta_ra_y, delta_dec_y)
    if raw_det > 0:
        cd = scale_deg_px * np.array([[-math.cos(theta), math.sin(theta)], [math.sin(theta), math.cos(theta)]])
    else:
        cd = scale_deg_px * np.array([[math.cos(theta), math.sin(theta)], [-math.sin(theta), math.cos(theta)]])

    header: dict[str, Any] = {
        "CTYPE1": "RA---TAN",
        "CTYPE2": "DEC--TAN",
        "CUNIT1": "deg",
        "CUNIT2": "deg",
        "CRVAL1": float(row["ctr_ra"]),
        "CRVAL2": float(row["ctr_dec"]),
        "CRPIX1": jpg_width_px / 2.0,
        "CRPIX2": jpg_height_px / 2.0,
        "CD1_1": float(cd[0, 0]),
        "CD1_2": float(cd[0, 1]),
        "CD2_1": float(cd[1, 0]),
        "CD2_2": float(cd[1, 1]),
        "EQUINOX": 2000.0,
        "RADESYS": "FK5",
    }
    if row.get("datetime"):
        header["DATE-OBS"] = row["datetime"]
    return header


def jpg_to_fits(jpg_path: Path, out_fits_path: Path, header: dict[str, Any] | None = None) -> Path:
    """Convert a plate JPG to a grayscale FITS for SExtractor, and save it.

    `header` (if given) is written into the FITS header as-is - e.g. a guess
    WCS, so SExtractor's own FITS_LDAC output carries it through for SCAMP.
    """
    with Image.open(jpg_path) as im:
        gray = np.array(im.convert("L"), dtype=np.float32)

    out_fits_path.parent.mkdir(parents=True, exist_ok=True)
    fits.PrimaryHDU(data=255.0 - gray, header=fits.Header(header or {})).writeto(out_fits_path, overwrite=True)
    return out_fits_path
