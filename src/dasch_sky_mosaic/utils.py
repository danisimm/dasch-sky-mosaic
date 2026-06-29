"""Shared utility functions used across dasch_sky_mosaic modules."""
from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, cast

import requests
from astropy.io import fits

LOG = logging.getLogger(__name__)


def _download_stream(url: str, dest: Path, overwrite: bool = True) -> int:
    """Stream a URL to disk. Returns file size in bytes."""
    if dest.exists() and not overwrite:
        return dest.stat().st_size
    dest.parent.mkdir(parents=True, exist_ok=True)
    with requests.get(url, stream=True, timeout=180.0) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                if chunk:
                    f.write(chunk)
    return dest.stat().st_size


def _unpacked_path(fz_path: Path) -> Path:
    """Return the .fits path that funpack would produce from a .fit.fz path."""
    out = fz_path.with_suffix("")
    if out.suffix.lower() == ".fit":
        return out.with_suffix(".fits")
    if out.suffix.lower() != ".fits":
        return out.with_name(out.name + ".fits")
    return out


def _decompress_fz(path: Path) -> Path:
    """Decompress a .fit.fz file to .fits via astropy, delete the original, and return the path."""
    if path.suffix.lower() != ".fz":
        return path

    out_path = _unpacked_path(path)

    if out_path.exists():
        path.unlink(missing_ok=True)
        return out_path

    LOG.info("Decompressing %s", path.name)
    written = False
    with fits.open(path, memmap=False) as hdul:
        for hdu in hdul:
            data = getattr(hdu, "data", None)
            if data is None or not hasattr(data, "ndim") or data.ndim < 2:
                continue
            primary = fits.PrimaryHDU(data=data.copy(), header=cast(Any, hdu).header.copy())
            fits.HDUList([primary]).writeto(out_path, overwrite=True)
            written = True
            break
    if not written:
        raise RuntimeError(f"no 2D image HDU found in {path.name}")
    path.unlink()
    return out_path
