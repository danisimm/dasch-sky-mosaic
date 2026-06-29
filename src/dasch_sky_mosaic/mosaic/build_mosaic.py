"""Build a stitched FITS mosaic from DASCH plates.

Run directly:
    python -m dasch_sky_mosaic.mosaic.build_mosaic --help

Or import build_mosaic() from dasch_sky_mosaic.mosaic.stitch for library use.
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
from pathlib import Path

from dotenv import load_dotenv

from dasch_sky_mosaic.call_sg import parse_cli_date_jd
from dasch_sky_mosaic.mosaic.discover import BuildConfig, Region
from dasch_sky_mosaic.mosaic.stitch import build_mosaic


def _parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build a stitched FITS mosaic from DASCH plates.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # Sky region
    p.add_argument("--ra", type=float, required=True, help="Region center RA (degrees)")
    p.add_argument("--dec", type=float, required=True, help="Region center Dec (degrees)")
    p.add_argument("--width", type=float, default=1.0, help="Region width (degrees)")
    p.add_argument("--height", type=float, default=None, help="Region height (degrees); defaults to --width")

    # Date filtering
    p.add_argument("--as-of-date", default=None, help="Only use plates observed on or before this date (ISO or JD)")
    p.add_argument("--earliest-date", default=None, help="Only use plates observed on or after this date (ISO or JD)")

    # Output
    p.add_argument("--output", required=True, help="Output FITS mosaic path")
    p.add_argument("--epoch-fits", default=None, help="Optional output path for per-pixel epoch FITS")
    p.add_argument("--manifest", default=None, help="Output manifest JSON path (defaults to output stem + _manifest.json)")
    p.add_argument("--session-dir", default="data/cache/dasch_session", help="Cache directory for downloaded files")
    p.add_argument("--overwrite", action="store_true", help="Overwrite existing output files")

    # Plate selection
    p.add_argument("--as-of-date-from-manifest", default=None, metavar="JSON", help="Load plate list from a prior manifest instead of running discovery")
    p.add_argument("--query-step", type=float, default=0.5, help="Grid spacing for plate discovery (degrees)")
    p.add_argument("--max-plates", type=int, default=None, help="Maximum number of plates to use")
    p.add_argument("--allow-multi-solution", action="store_true", help="Allow plates with multiple WCS solutions")

    # Image parameters
    p.add_argument("--pixel-scale", type=float, default=None, help="Output pixel scale (arcsec/pixel); defaults to native plate scale")
    p.add_argument("--projection", choices=["TAN", "CAR"], default="TAN", help="WCS projection type")
    p.add_argument("--binning", type=int, choices=[1, 16], default=16, help="DASCH mosaic binning level")

    # Background matching
    p.add_argument("--subtract-background", action="store_true", help="Apply per-plate background subtraction and global matching")

    # API
    p.add_argument("--api-base", default="auto", help="Starglass API base URL (auto/public/full or explicit URL)")
    p.add_argument("--api-key", default=None, help="Starglass API key (or set STARGLASS_API_KEY env var)")

    # Logging
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    load_dotenv()
    args = _parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(levelname)s %(message)s",
        stream=sys.stderr,
    )

    api_key = args.api_key or os.environ.get("STARGLASS_API_KEY")
    output = Path(args.output)
    manifest = Path(args.manifest) if args.manifest else output.with_name(output.stem + "_manifest.json")
    session_root = Path(args.session_dir)
    height = args.height if args.height is not None else args.width

    config = BuildConfig(
        region=Region(ra_deg=args.ra, dec_deg=args.dec, width_deg=args.width, height_deg=height),
        as_of_jd=parse_cli_date_jd(args.as_of_date),
        earliest_jd=parse_cli_date_jd(args.earliest_date),
        session_root=session_root,
        output_fits=output,
        epoch_fits=Path(args.epoch_fits) if args.epoch_fits else None,
        manifest_json=manifest,
        pixel_scale_arcsec=args.pixel_scale,
        projection=args.projection,
        binning=args.binning,
        query_step_deg=args.query_step,
        api_base=args.api_base,
        api_key=api_key,
        subtract_background=args.subtract_background,
        allow_multi_solution_plates=args.allow_multi_solution,
        delete_base_mosaics=False,
        overwrite=args.overwrite,
        max_plates=args.max_plates,
        from_manifest=Path(args.as_of_date_from_manifest) if args.as_of_date_from_manifest else None,
    )

    result = build_mosaic(config)
    print(f"Done. {result['plates'].__len__()} plates → {output}")


if __name__ == "__main__":
    main()
