from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

from dotenv import load_dotenv

from dasch_sky_mosaic.pipeline import BuildConfig, Region, build_mosaic, parse_cli_date_jd


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dasch-sky-mosaic",
        description="Build a time-filtered, WCS-preserving mosaic from DASCH full-plate mosaics.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build", help="Query DASCH and build a stitched FITS mosaic")
    build_parser.add_argument("--ra-deg", type=float, required=True, help="Mosaic center right ascension in degrees")
    build_parser.add_argument("--dec-deg", type=float, required=True, help="Mosaic center declination in degrees")
    build_parser.add_argument("--width-deg", type=float, required=True, help="Width of the target sky region in degrees")
    build_parser.add_argument("--height-deg", type=float, required=True, help="Height of the target sky region in degrees")
    build_parser.add_argument("--as-of-date", help="Only include plates observed on or before this date (ISO, e.g. 1954-03-01)")
    build_parser.add_argument("--earliest-date", help="Only include plates observed on or after this date (ISO, e.g. 1890-01-01)")
    build_parser.add_argument(
        "--pixel-scale-arcsec",
        type=float,
        help="Output pixel scale in arcseconds per pixel (default: native plate scale of selected mosaics)",
    )
    build_parser.add_argument("--projection", choices=["TAN", "CAR"], default="TAN", help="Output WCS projection")
    build_parser.add_argument("--binning", choices=[1, 16], type=int, default=16, help="DASCH mosaic binning level")
    build_parser.add_argument("--query-step-deg", type=float, default=5.0, help="Sampling step for plate discovery across the requested region")
    build_parser.add_argument("--session-root", type=Path, default=Path("data/cache/dasch_session"), help="Local daschlab cache directory")
    build_parser.add_argument("--output-fits", type=Path, required=True, help="Path for the stitched output FITS")
    build_parser.add_argument("--epoch-fits", type=Path, help="Optional path for a per-pixel observation-date FITS (JD values)")
    build_parser.add_argument("--manifest-json", type=Path, required=True, help="Path for the JSON build manifest")
    build_parser.add_argument("--api-base", default="auto", help="API base alias: auto, public, full, or a full custom base URL")
    build_parser.add_argument("--api-key-env", default="DASCHLAB_API_KEY", help="Environment variable containing the Starglass/DASCHLAB API key")
    build_parser.add_argument("--allow-multi-solution-plates", action="store_true", help="Include plates with multiple WCS solutions")
    build_parser.add_argument("--preserve-native-background", action="store_true", help="Disable per-plate median background subtraction")
    build_parser.add_argument("--delete-base-mosaics", action="store_true", help="Delete cached base mosaics after building to free disk space")
    build_parser.add_argument("--max-plates", type=int, help="Optional hard cap on the number of selected plates")
    build_parser.add_argument("--overwrite", action="store_true", help="Allow overwriting existing output files")
    build_parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    load_dotenv()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(message)s", force=True)
    log = logging.getLogger(__name__)

    config = BuildConfig(
        region=Region(
            ra_deg=args.ra_deg,
            dec_deg=args.dec_deg,
            width_deg=args.width_deg,
            height_deg=args.height_deg,
        ),
        as_of_jd=parse_cli_date_jd(args.as_of_date),
        earliest_jd=parse_cli_date_jd(args.earliest_date),
        session_root=args.session_root,
        output_fits=args.output_fits,
        epoch_fits=args.epoch_fits,
        manifest_json=args.manifest_json,
        pixel_scale_arcsec=args.pixel_scale_arcsec,
        projection=args.projection,
        binning=args.binning,
        query_step_deg=args.query_step_deg,
        api_base=args.api_base,
        api_key=os.getenv(args.api_key_env),
        subtract_background=not args.preserve_native_background,
        allow_multi_solution_plates=args.allow_multi_solution_plates,
        delete_base_mosaics=args.delete_base_mosaics,
        overwrite=args.overwrite,
        max_plates=args.max_plates,
    )
    log.info(
        "Starting DASCH mosaic build: center=(%.6f, %.6f) size=%.3fx%.3f deg binning=%d",
        args.ra_deg,
        args.dec_deg,
        args.width_deg,
        args.height_deg,
        args.binning,
    )
    manifest = build_mosaic(config)
    log.info(
        "Build complete: %d plate(s) -> %s",
        len(manifest.get("plates", [])),
        args.output_fits,
    )
    return 0