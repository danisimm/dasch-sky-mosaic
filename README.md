# DASCH Sky Mosaic Builder

This workspace builds WCS-preserving sky mosaics from DASCH full-plate images.
It uses DASCH/Starglass APIs to discover candidate plates and `daschlab` to download the standardized value-added mosaic FITS files, because the raw mosaic APIs do not provide the best available WCS metadata.

## What the pipeline does

- Queries the DASCH exposure catalog over a user-defined sky region by sampling a grid of sky positions.
- Uses only `queryexps` response fields (`obsDate`, `solnum`, `wcssource`) to pick candidates, avoiding extra per-plate metadata calls.
- At each sampled sky position, selects the single most recent eligible plate.
- Downloads DASCH value-added mosaics through `daschlab.Session.mosaic()`.
- Reprojects each selected plate onto a user-defined output WCS.
- Stitches with a last-write-wins strategy (oldest to newest), so each output pixel reflects the most recent plate covering that point.
- Writes a JSON manifest describing the selected inputs and observation dates.

## Important scientific limitations

- DASCH exposure discovery is point-based. For an extended region, this project queries a configurable grid and unions the per-point winners. Smaller grid spacing improves completeness at the cost of more API calls.
- By default, only plates with exactly one astrometric solution are included. This avoids misprojecting multi-exposure plates, where one plate image can contain multiple distinct sky solutions.
- The default combination mode is "most recent wins" after optional per-plate median background subtraction. This is suitable for visualization and chronology tracking, not precision photometry.
- The output FITS keeps accurate WCS in the target projection, but this project does not yet generate WWT TOAST tile pyramids. The output mosaic is intended to be the scientifically valid intermediate product for later tiling.

## Installation

```powershell
python -m venv .venv
.venv\Scripts\Activate.ps1
python -m pip install --upgrade pip
python -m pip install -e .
```

If you have a Starglass API key, set it to raise API rate limits:

```powershell
$env:DASCHLAB_API_KEY = "your-api-key"
```

Alternatively, place the key in a `.env` file at the project root (gitignored):

```
DASCHLAB_API_KEY=your-api-key
```

The CLI loads `.env` automatically on startup. Without an API key the code uses the public (rate-limited) API.

## Example

This example builds a 6-by-6 degree TAN-projected mosaic around M31 using all eligible plates from 1910 up to 1930, with an epoch map:

```powershell
python -m dasch_sky_mosaic build `
  --ra-deg 10.6847083 `
  --dec-deg 41.26875 `
  --width-deg 6 `
  --height-deg 6 `
  --earliest-date 1910-01-01 `
  --as-of-date 1930-12-31 `
  --pixel-scale-arcsec 60 `
  --output-fits data/output/m31_1910_1930.fits `
  --epoch-fits data/output/m31_1910_1930_epoch.fits `
  --delete-base-mosaics `
  --manifest-json data/output/m31_1910_1930_manifest.json
```

## Key options

- `--binning {1,16}`: DASCH mosaic binning level. `16` is dramatically smaller and should be the default for large-area work.
- `--query-step-deg`: spacing of the discovery grid used to find plates intersecting the requested sky region.
- `--projection {TAN,CAR}`: output FITS projection.
- `--as-of-date` and `--earliest-date`: optional date bounds for candidate exposures. Omit both to build the newest-available mosaic.
- `--allow-multi-solution-plates`: include plates with multiple WCS solutions. This can produce incorrect reprojections and is disabled by default.
- `--preserve-native-background`: skip median background subtraction before coaddition.
- `--delete-base-mosaics`: remove cached base mosaics after build to free disk space.

Current implementation note:
- FITS mode is fully supported.

## Output products

- Main mosaic FITS with the requested output WCS.
- Optional epoch FITS (`--epoch-fits`) with per-pixel Julian Date values for the winning plate.
- JSON manifest recording selected plates, most-recent observation dates, and local cache paths.

## Next step for Worldwide Telescope

Once the stitched FITS looks right, convert it into a TOAST tile pyramid using a WWT-compatible tiler such as `toasty` or `wwt_data_formats`. This script is focused on the plate-selection, download, reprojection, and coaddition stage where WCS integrity matters most.