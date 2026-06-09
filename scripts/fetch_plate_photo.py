#!/usr/bin/env python3
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from urllib.parse import urlparse

import requests

from dasch_sky_mosaic.fetch import StarglassClient, jd_to_iso


def filename_from_url(url: str, fallback_stem: str) -> str:
    path = urlparse(url).path
    leaf = Path(path).name
    if not leaf:
        return f"{fallback_stem}.jpg"
    return leaf


def download(url: str, dest: Path, overwrite: bool = False) -> int:
    dest.parent.mkdir(parents=True, exist_ok=True)
    if dest.exists() and not overwrite:
        return dest.stat().st_size
    with requests.get(url, stream=True, timeout=180.0) as r:
        r.raise_for_status()
        with dest.open("wb") as f:
            for chunk in r.iter_content(chunk_size=65536):
                if chunk:
                    f.write(chunk)
    return dest.stat().st_size


def main(argv: list[str] | None = None) -> int:
    argv = argv or sys.argv[1:]
    if not argv:
        print("Usage: fetch_plate_photo.py <plate_id>")
        return 2
    plate_id = argv[0].strip().lower()
    api_key = os.getenv("DASCHLAB_API_KEY")
    client = StarglassClient(api_base="auto", api_key=api_key)

    payload = client.get_plate(plate_id)
    images = [
        c for c in (payload.get("plate_images") or [])
        if not bool(c.get("thumbnail", False)) and str(c.get("portion", "")).lower() == "all"
    ]
    if not images:
        print(f"No full-plate images found for {plate_id}")
        return 1
    best = max(images, key=lambda e: int(e.get("thumbnail_ratio") or 0))
    url = str(best.get("url") or "").strip()
    if not url:
        print(f"No URL available for {plate_id}")
        return 1

    out_dir = Path("data/cache/dasch_session/plate_photos")
    filename = filename_from_url(url, fallback_stem=f"{plate_id}_plate_all")
    dest = out_dir / filename
    try:
        nbytes = download(url, dest, overwrite=True)
    except Exception as exc:
        print(f"Download failed: {exc}")
        return 1

    obs_jd = None
    for k in ("obsDate", "obs_date", "expdate", "expDate"):
        if payload.get(k):
            try:
                obs_jd = float(payload[k])
            except Exception:
                pass
    manifest = {
        "workflow": "plate-photo-download",
        "region": None,
        "as_of_date": None,
        "earliest_date": None,
        "query_step_deg": None,
        "n_candidates": 1,
        "n_downloaded": 1,
        "photos": [
            {
                "plate_id": plate_id,
                "obs_date": jd_to_iso(obs_jd) if obs_jd else None,
                "obs_date_jd": obs_jd,
                "status": "downloaded",
                "portion": best.get("portion"),
                "image_type": best.get("image_type"),
                "thumbnail": bool(best.get("thumbnail", False)),
                "local_photo_path": str(dest),
                "bytes": nbytes,
            }
        ],
    }

    out_json = Path("data/output") / f"{plate_id}_photos.json"
    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(out_json)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
