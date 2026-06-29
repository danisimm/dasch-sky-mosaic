"""Shared utility functions used across dasch_sky_mosaic modules."""
from __future__ import annotations

from pathlib import Path

import requests


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
