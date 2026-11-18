# dasch-sky-mosaic

Personal research tools for working with DASCH photographic plate data.
Two distinct workflows: a WCS batch pipeline for generating WTML placement data,
and an all-sky mosaic pipeline for stitching plates into large FITS images.

## Installation

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
pip install -e .
```

If you have a Starglass API key, set it to use the authenticated (higher rate-limit) API:

```powershell
$env:STARGLASS_API_KEY = "your-api-key"
```

Or place it in a `.env` file at the project root (gitignored):

```
STARGLASS_API_KEY=your-api-key
```

Credentials and API key files live in `credentials/` (gitignored).

---

## Workflow 1: WCS batch solve

Given a list of plate IDs, download the plate JPG and its binning-16 FITS mosaic,
solve WCS for the JPG by aligning it to the FITS (astroalign star-pattern matching),
and write per-plate JSON for WTML generation. Intended to run at scale across 400k+ plates.

```powershell
python -m dasch_sky_mosaic.pipeline a11740 br00130
```

Settings (download method `s3`/`sg`, cache and output directories) are set in the
configuration block at the bottom of `pipeline.py`.

### Output

Each plate produces `data/output/wcs/{plate_id}_wcs.json` containing:
- `wcs_header`: FITS WCS keywords for the plate JPG
- `image_width_px`, `image_height_px`
- `alignment_meta`: solver diagnostics (backend, parity, matches, transform scale/rotation)

`wtml_gen.py` converts these records into WWT placement fields and WTML. In production,
WTML generation is meant to happen server-side in an AWS Lambda function.

A JPG-only SExtractor + SCAMP solver was tried and shelved; see the
`archive/scamp-single-image` branch (`SCAMP_ARCHIVE.md`).

---

## Workflow 2: All-sky mosaic

Discover plates covering a sky region and time range, download their FITS mosaics,
reproject onto a common WCS, and stitch into a single output FITS.

```powershell
python -m dasch_sky_mosaic.mosaic.build_mosaic `
    --ra 83.8 --dec -5.4 `
    --width 2.0 `
    --as-of-date 1950-01-01 `
    --output data/output/orion_1950.fits
```

### Key options

- `--ra`, `--dec`: region center (degrees, required).
- `--width`, `--height`: region size in degrees (height defaults to width).
- `--as-of-date`, `--earliest-date`: ISO or JD date bounds for candidate plates.
- `--output`: output FITS mosaic path (required).
- `--epoch-fits`: optional output for per-pixel Julian Date values.
- `--manifest`: output manifest JSON (defaults to `{output stem}_manifest.json`).
- `--session-dir`: cache directory for downloaded files (default: `data/cache/dasch_session`).
- `--pixel-scale`: output pixel scale in arcsec/pixel (defaults to native plate scale).
- `--projection {TAN,CAR}`: WCS projection (default: `TAN`).
- `--binning {1,16}`: DASCH mosaic binning level (default: `16`; much smaller downloads).
- `--query-step`: grid spacing for plate discovery in degrees (default: `0.5`).
- `--max-plates`: cap the number of plates used.
- `--allow-multi-solution`: include plates with multiple WCS solutions (off by default).
- `--subtract-background`: apply per-plate background subtraction and global matching.
- `--overwrite`: overwrite existing output files.
- `--api-key`: Starglass API key (or set `STARGLASS_API_KEY` env var).

### Output

- Main mosaic FITS with the requested WCS.
- Optional epoch FITS with per-pixel Julian Date of the winning plate.
- JSON manifest recording selected plates and local cache paths.

---

## Module structure

```
src/dasch_sky_mosaic/
  call_sg.py         # Starglass API client; JPG and FITS download via API
  call_s3.py         # Direct S3 downloads (JPG and binning-16 FITS)
  pipeline.py        # Runnable batch script: download + solve + save JSON
  wcs_solve.py       # Astroalign WCS solver (JPG aligned to FITS mosaic)
  wtml_gen.py        # Pure WTML XML generation from WCS placement dicts
  mosaic/
    discover.py      # Region/BuildConfig dataclasses; plate grid discovery via queryexps
    background.py    # Montage-style global background matching
    plate_photos.py  # Region-based plate photo (JPG) discovery and download
    stitch.py        # Reprojection and coaddition pipeline
    build_mosaic.py  # Runnable mosaic builder
```
