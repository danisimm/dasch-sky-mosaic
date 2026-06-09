# This file was generated with the assistance of GitHub Copilot (AI).
from __future__ import annotations

import argparse
import logging
import os
from pathlib import Path

from dotenv import load_dotenv

from dasch_sky_mosaic.fetch import BuildConfig, Region, parse_cli_date_jd
from dasch_sky_mosaic.pipeline import build_mosaic
from dasch_sky_mosaic.plate_photos import PlatePhotoConfig, discover_and_download_plate_photos
from dasch_sky_mosaic.wtml import WtmlBuildConfig, build_wtml


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
    build_parser.add_argument(
        "--from-manifest",
        type=Path,
        help="Path to an existing build manifest JSON; bypasses API discovery and download, using the cached plate files listed in the manifest",
    )
    build_parser.add_argument("--overwrite", action="store_true", help="Allow overwriting existing output files")
    build_parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")

    photos_parser = subparsers.add_parser(
        "plate-photos",
        help="Discover candidate plates for a sky region and download source JPG plate photos",
    )
    photos_parser.add_argument("--ra-deg", type=float, required=True, help="Region center right ascension in degrees")
    photos_parser.add_argument("--dec-deg", type=float, required=True, help="Region center declination in degrees")
    photos_parser.add_argument("--width-deg", type=float, required=True, help="Width of the target sky region in degrees")
    photos_parser.add_argument("--height-deg", type=float, required=True, help="Height of the target sky region in degrees")
    photos_parser.add_argument("--as-of-date", help="Only include plates observed on or before this date (ISO)")
    photos_parser.add_argument("--earliest-date", help="Only include plates observed on or after this date (ISO)")
    photos_parser.add_argument("--query-step-deg", type=float, default=5.0, help="Sampling step for plate discovery across the requested region")
    photos_parser.add_argument("--api-base", default="auto", help="API base alias: auto, public, full, or a full custom base URL")
    photos_parser.add_argument("--api-key-env", default="DASCHLAB_API_KEY", help="Environment variable containing the Starglass/DASCHLAB API key")
    photos_parser.add_argument("--allow-multi-solution-plates", action="store_true", help="Include plates with multiple WCS solutions")
    photos_parser.add_argument("--max-plates", type=int, help="Optional hard cap on the number of selected plates")
    photos_parser.add_argument("--output-dir", type=Path, default=Path("data/cache/dasch_session/plate_photos"), help="Directory to store downloaded plate photo JPG files")
    photos_parser.add_argument("--manifest-json", type=Path, required=True, help="Path for the output JSON manifest describing downloaded photos")
    photos_parser.add_argument("--overwrite", action="store_true", help="Allow overwriting existing photo files/manifest")
    photos_parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")

    wtml_parser = subparsers.add_parser(
        "wtml",
        help="Photo+WCS->WTML pipeline (supports fresh discovery or existing photo manifest)",
    )
    wtml_parser.add_argument("--output-wtml", type=Path, required=True, help="Output WTML path")
    wtml_parser.add_argument("--output-json", type=Path, required=True, help="Output pairing/report JSON path")
    wtml_parser.add_argument("--photo-manifest-json", type=Path, help="Optional existing plate-photos manifest; if omitted, discovery runs")
    wtml_parser.add_argument("--ra-deg", type=float, help="Region center right ascension in degrees")
    wtml_parser.add_argument("--dec-deg", type=float, help="Region center declination in degrees")
    wtml_parser.add_argument("--width-deg", type=float, help="Width of the target sky region in degrees")
    wtml_parser.add_argument("--height-deg", type=float, help="Height of the target sky region in degrees")
    wtml_parser.add_argument("--as-of-date", help="Only include plates observed on or before this date (ISO)")
    wtml_parser.add_argument("--earliest-date", help="Only include plates observed on or after this date (ISO)")
    wtml_parser.add_argument("--query-step-deg", type=float, default=5.0, help="Sampling step for discovery when photo manifest is omitted")
    wtml_parser.add_argument("--api-base", default="auto", help="API base alias: auto, public, full, or custom URL")
    wtml_parser.add_argument("--api-key-env", default="DASCHLAB_API_KEY", help="Environment variable containing API key")
    wtml_parser.add_argument("--allow-multi-solution-plates", action="store_true", help="Include plates with multiple WCS solutions")
    wtml_parser.add_argument("--max-plates", type=int, help="Optional hard cap on discovered plates")
    wtml_parser.add_argument("--session-root", type=Path, default=Path("data/cache/dasch_session"), help="daschlab cache root for FITS WCS retrieval")
    wtml_parser.add_argument("--photo-output-dir", type=Path, default=Path("data/cache/dasch_session/plate_photos"), help="Photo download directory when discovery runs")
    wtml_parser.add_argument("--overwrite", action="store_true", help="Allow overwriting output artifacts")
    wtml_parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    load_dotenv()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(message)s", force=True)
    log = logging.getLogger(__name__)

    if args.command == "build":
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
            from_manifest=args.from_manifest,
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
    elif args.command == "plate-photos":
        config = PlatePhotoConfig(
            region=Region(
                ra_deg=args.ra_deg,
                dec_deg=args.dec_deg,
                width_deg=args.width_deg,
                height_deg=args.height_deg,
            ),
            as_of_jd=parse_cli_date_jd(args.as_of_date),
            earliest_jd=parse_cli_date_jd(args.earliest_date),
            query_step_deg=args.query_step_deg,
            api_base=args.api_base,
            api_key=os.getenv(args.api_key_env),
            allow_multi_solution_plates=args.allow_multi_solution_plates,
            max_plates=args.max_plates,
            output_dir=args.output_dir,
            manifest_json=args.manifest_json,
            overwrite=args.overwrite,
        )
        log.info(
            "Starting DASCH plate-photo download: center=(%.6f, %.6f) size=%.3fx%.3f deg",
            args.ra_deg,
            args.dec_deg,
            args.width_deg,
            args.height_deg,
        )
        manifest = discover_and_download_plate_photos(config)
        log.info(
            "Photo download complete: %d/%d plate(s) have local photos -> %s",
            manifest.get("n_downloaded", 0),
            manifest.get("n_candidates", 0),
            args.manifest_json,
        )
    elif args.command == "wtml":
        region = None
        if args.photo_manifest_json is None:
            required = [args.ra_deg, args.dec_deg, args.width_deg, args.height_deg]
            if any(v is None for v in required):
                parser.error("wtml requires either --photo-manifest-json or full region args (--ra-deg/--dec-deg/--width-deg/--height-deg)")
            region = Region(
                ra_deg=float(args.ra_deg),
                dec_deg=float(args.dec_deg),
                width_deg=float(args.width_deg),
                height_deg=float(args.height_deg),
            )

        config = WtmlBuildConfig(
            output_wtml=args.output_wtml,
            output_json=args.output_json,
            photo_manifest_json=args.photo_manifest_json,
            region=region,
            as_of_date=args.as_of_date,
            earliest_date=args.earliest_date,
            query_step_deg=args.query_step_deg,
            api_base=args.api_base,
            api_key=os.getenv(args.api_key_env),
            allow_multi_solution_plates=args.allow_multi_solution_plates,
            max_plates=args.max_plates,
            session_root=args.session_root,
            photo_output_dir=args.photo_output_dir,
            overwrite=args.overwrite,
        )
        log.info("Starting WTML pipeline")
        result = build_wtml(config)
        log.info(
            "WTML complete: paired=%d rejected=%d -> %s",
            result.get("n_paired", 0),
            result.get("n_rejected", 0),
            args.output_wtml,
        )
    else:  # pragma: no cover
        parser.error(f"unsupported command: {args.command}")
    return 0