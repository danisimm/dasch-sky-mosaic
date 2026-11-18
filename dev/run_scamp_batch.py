"""Run the SExtractor+SCAMP WCS pipeline on every plate already cached in data/cache/photos.

Run from WSL (repo root):
    python dev/run_scamp_batch.py
"""
import subprocess
import sys
from pathlib import Path

PHOTO_DIR = Path("data/cache/photos")
SUFFIX = "_pphoto_all.jpg"

plate_ids = sorted(p.name[: -len(SUFFIX)] for p in PHOTO_DIR.glob(f"*{SUFFIX}"))
print(f"Found {len(plate_ids)} cached plates: {plate_ids}")

subprocess.run([sys.executable, "-m", "dasch_sky_mosaic.wcs_solve_scamp", *plate_ids], check=False)
