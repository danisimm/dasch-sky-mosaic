# This file was generated with the assistance of GitHub Copilot (AI).
from __future__ import annotations

import argparse
import logging
import os
import re
from pathlib import Path

from astropy.coordinates import SkyCoord
from dotenv import load_dotenv

from dasch_sky_mosaic.fetch import BuildConfig, Region, parse_cli_date_jd
from dasch_sky_mosaic.pipeline import build_mosaic
from dasch_sky_mosaic.plate_photos import PlatePhotoConfig, discover_and_download_plate_photos
from dasch_sky_mosaic.wtml import WtmlBuildConfig, WtmlPlateSolveConfig, build_wtml, build_wtml_from_plate_solve


def _resolve_target(target_name: str) -> tuple[float, float]:
    """Resolve a target name to RA/Dec using SkyCoord.from_name (Simbad)."""
    try:
        coord = SkyCoord.from_name(target_name)
        return float(coord.ra.deg), float(coord.dec.deg)
    except Exception as e:
        raise ValueError(f"Could not resolve target '{target_name}' via Simbad: {e}") from e


def _make_output_slug(target_name: str | None, ra_deg: float, dec_deg: float) -> str:
    """Generate a slug for auto-naming output files."""
    if target_name:
        # Sanitize: lowercase, replace spaces with underscores, remove non-alphanumeric except underscore
        slug = re.sub(r'[^a-z0-9_]', '', target_name.lower().replace(' ', '_'))
        return slug
    else:
        # Use coordinates: "10.68_41.27" format (no sign prefix)
        return f"{ra_deg:.2f}_{abs(dec_deg):.2f}"


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="dasch-sky-mosaic",
        description="Build a time-filtered, WCS-preserving mosaic from DASCH full-plate mosaics.",
    )
    subparsers = parser.add_subparsers(dest="command", required=True)

    build_parser = subparsers.add_parser("build", help="Query DASCH and build a stitched FITS mosaic")
    build_parser.add_argument("--target", help="Target name (e.g. 'Polaris Australis'); resolves via Simbad")
    build_parser.add_argument("--ra-deg", type=float, help="Mosaic center right ascension in degrees")
    build_parser.add_argument("--dec-deg", type=float, help="Mosaic center declination in degrees")
    build_parser.add_argument("--size-deg", type=float, default=1.0, help="Size of the region in degrees (sets both width and height)")
    build_parser.add_argument("--width-deg", type=float, help="Width of the target sky region in degrees (overrides --size-deg)")
    build_parser.add_argument("--height-deg", type=float, help="Height of the target sky region in degrees (overrides --size-deg)")
    build_parser.add_argument("--as-of-date", help="Only include plates observed on or before this date (ISO, e.g. 1954-03-01)")
    build_parser.add_argument("--earliest-date", help="Only include plates observed on or after this date (ISO, e.g. 1890-01-01)")
    build_parser.add_argument(
        "--pixel-scale-arcsec",
        type=float,
        help="Output pixel scale in arcseconds per pixel (default: native plate scale of selected mosaics)",
    )
    build_parser.add_argument("--projection", choices=["TAN", "CAR"], default="TAN", help="Output WCS projection")
    build_parser.add_argument("--binning", choices=[1, 16], type=int, default=16, help="DASCH mosaic binning level")
    build_parser.add_argument("--query-step-deg", type=float, default=999.0, help="Sampling step for plate discovery (default 999.0 queries only center point)")
    build_parser.add_argument("--session-root", type=Path, default=Path("data/cache/dasch_session"), help="Local daschlab cache directory")
    build_parser.add_argument("--output-fits", type=Path, help="Path for the stitched output FITS (auto-generated if omitted)")
    build_parser.add_argument("--epoch-fits", type=Path, help="Optional path for a per-pixel observation-date FITS (JD values)")
    build_parser.add_argument("--manifest-json", type=Path, help="Path for the JSON build manifest (auto-generated if omitted)")
    build_parser.add_argument("--api-base", default="auto", help="API base alias: auto, public, full, or a full custom base URL")
    build_parser.add_argument("--api-key-env", default="DASCHLAB_API_KEY", help="Environment variable containing the Starglass/DASCHLAB API key")
    build_parser.add_argument("--allow-multi-solution-plates", action="store_true", help="Include plates with multiple WCS solutions")
    build_parser.add_argument("--preserve-native-background", action="store_true", help="Disable per-plate median background subtraction")
    build_parser.add_argument("--delete-base-mosaics", action="store_true", help="Delete cached base mosaics after building to free disk space")
    build_parser.add_argument("--max-plates", type=int, default=1, help="Hard cap on the number of selected plates (default 1)")
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
    photos_parser.add_argument("--target", help="Target name (e.g. 'Polaris Australis'); resolves via Simbad")
    photos_parser.add_argument("--ra-deg", type=float, help="Region center right ascension in degrees")
    photos_parser.add_argument("--dec-deg", type=float, help="Region center declination in degrees")
    photos_parser.add_argument("--size-deg", type=float, default=1.0, help="Size of the region in degrees (sets both width and height)")
    photos_parser.add_argument("--width-deg", type=float, help="Width of the target sky region in degrees (overrides --size-deg)")
    photos_parser.add_argument("--height-deg", type=float, help="Height of the target sky region in degrees (overrides --size-deg)")
    photos_parser.add_argument("--as-of-date", help="Only include plates observed on or before this date (ISO)")
    photos_parser.add_argument("--earliest-date", help="Only include plates observed on or after this date (ISO)")
    photos_parser.add_argument("--query-step-deg", type=float, default=999.0, help="Sampling step for plate discovery (default 999.0 queries only center point)")
    photos_parser.add_argument("--api-base", default="auto", help="API base alias: auto, public, full, or a full custom base URL")
    photos_parser.add_argument("--api-key-env", default="DASCHLAB_API_KEY", help="Environment variable containing the Starglass/DASCHLAB API key")
    photos_parser.add_argument("--allow-multi-solution-plates", action="store_true", help="Include plates with multiple WCS solutions")
    photos_parser.add_argument("--max-plates", type=int, default=1, help="Hard cap on the number of selected plates (default 1)")
    photos_parser.add_argument("--output-dir", type=Path, default=Path("data/cache/dasch_session/plate_photos"), help="Directory to store downloaded plate photo JPG files")
    photos_parser.add_argument("--manifest-json", type=Path, help="Path for the output JSON manifest (auto-generated if omitted)")
    photos_parser.add_argument("--overwrite", action="store_true", help="Allow overwriting existing photo files/manifest")
    photos_parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")

    wtml_parser = subparsers.add_parser(
        "wtml",
        help="Photo+WCS->WTML pipeline (supports fresh discovery or existing photo manifest)",
    )
    wtml_parser.add_argument("--target", help="Target name (e.g. 'Polaris Australis'); resolves via Simbad")
    wtml_parser.add_argument("--output-wtml", type=Path, help="Output WTML path (auto-generated if omitted)")
    wtml_parser.add_argument("--output-json", type=Path, help="Output pairing/report JSON path (auto-generated if omitted)")
    wtml_parser.add_argument("--photo-manifest-json", type=Path, help="Optional existing plate-photos manifest; if omitted, discovery runs")
    wtml_parser.add_argument("--ra-deg", type=float, help="Region center right ascension in degrees")
    wtml_parser.add_argument("--dec-deg", type=float, help="Region center declination in degrees")
    wtml_parser.add_argument("--size-deg", type=float, default=1.0, help="Size of the region in degrees (sets both width and height)")
    wtml_parser.add_argument("--width-deg", type=float, help="Width of the target sky region in degrees (overrides --size-deg)")
    wtml_parser.add_argument("--height-deg", type=float, help="Height of the target sky region in degrees (overrides --size-deg)")
    wtml_parser.add_argument("--as-of-date", help="Only include plates observed on or before this date (ISO)")
    wtml_parser.add_argument("--earliest-date", help="Only include plates observed on or after this date (ISO)")
    wtml_parser.add_argument("--query-step-deg", type=float, default=999.0, help="Sampling step for discovery (default 999.0 queries only center point)")
    wtml_parser.add_argument("--api-base", default="auto", help="API base alias: auto, public, full, or custom URL")
    wtml_parser.add_argument("--api-key-env", default="DASCHLAB_API_KEY", help="Environment variable containing API key")
    wtml_parser.add_argument("--allow-multi-solution-plates", action="store_true", help="Include plates with multiple WCS solutions")
    wtml_parser.add_argument("--max-plates", type=int, default=1, help="Hard cap on discovered plates (default 1)")
    wtml_parser.add_argument("--session-root", type=Path, default=Path("data/cache/dasch_session"), help="daschlab cache root for FITS WCS retrieval")
    wtml_parser.add_argument("--photo-output-dir", type=Path, default=Path("data/cache/dasch_session/plate_photos"), help="Photo download directory when discovery runs")
    wtml_parser.add_argument("--s3-region", default="us-east-1", help="AWS region for S3 presigned URLs (default us-east-1)")
    wtml_parser.add_argument("--overwrite", action="store_true", help="Allow overwriting output artifacts")
    wtml_parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")

    wtml_solve_parser = subparsers.add_parser(
        "wtml-solve-plate",
        help="Build WTML for a single plate_id using the astrometry.net API on the full-plate photo",
    )
    wtml_solve_parser.add_argument("--plate-id", required=True, help="DASCH plate identifier (e.g. dnb06613)")
    wtml_solve_parser.add_argument("--output-wtml", type=Path, help="Output WTML path (auto-generated if omitted)")
    wtml_solve_parser.add_argument("--output-json", type=Path, help="Output report JSON path (auto-generated if omitted)")
    wtml_solve_parser.add_argument("--photo-output-dir", type=Path, default=Path("data/cache/dasch_session/plate_photos"), help="Directory for downloaded plate photos")
    wtml_solve_parser.add_argument("--solver-work-dir", type=Path, default=Path("data/cache/dasch_session/astrometry"), help="Working directory for astrometry.net outputs")
    wtml_solve_parser.add_argument("--astrometry-api-key", help="Astrometry.net API key (if omitted, read from --astrometry-api-key-env)")
    wtml_solve_parser.add_argument("--astrometry-api-key-env", default="ASTROMETRY_NET_API_KEY", help="Environment variable containing astrometry.net API key")
    wtml_solve_parser.add_argument("--astrometry-api-base", default="https://nova.astrometry.net/api", help="Astrometry.net API base URL")
    wtml_solve_parser.add_argument("--solver-timeout-sec", type=float, default=600.0, help="End-to-end solve timeout in seconds")
    wtml_solve_parser.add_argument("--solver-upload-timeout-sec", type=float, default=120.0, help="HTTP timeout for API calls in seconds")
    wtml_solve_parser.add_argument("--solver-poll-interval-sec", type=float, default=5.0, help="Polling interval for submission/job status checks")
    wtml_solve_parser.add_argument("--ra-hint-deg", type=float, help="Optional RA hint (deg)")
    wtml_solve_parser.add_argument("--dec-hint-deg", type=float, help="Optional Dec hint (deg)")
    wtml_solve_parser.add_argument("--radius-hint-deg", type=float, help="Optional search radius hint (deg)")
    wtml_solve_parser.add_argument("--scale-low-arcsec-per-pix", type=float, help="Optional lower pixel scale hint")
    wtml_solve_parser.add_argument("--scale-high-arcsec-per-pix", type=float, help="Optional upper pixel scale hint")
    wtml_solve_parser.add_argument("--api-base", default="auto", help="API base alias: auto, public, full, or custom URL")
    wtml_solve_parser.add_argument("--api-key-env", default="DASCHLAB_API_KEY", help="Environment variable containing API key")
    wtml_solve_parser.add_argument("--s3-region", default="us-east-1", help="AWS region for S3 presigned URLs (default us-east-1)")
    wtml_solve_parser.add_argument("--overwrite", action="store_true", help="Allow overwriting output artifacts")
    wtml_solve_parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"], help="Logging verbosity")
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)
    load_dotenv()
    logging.basicConfig(level=getattr(logging, args.log_level), format="%(levelname)s %(message)s", force=True)
    log = logging.getLogger(__name__)

    if args.command == "build":
        # Resolve target name if provided
        if args.target and (args.ra_deg or args.dec_deg):
            log.error("--target and --ra-deg/--dec-deg are mutually exclusive")
            return 1
        
        if args.target:
            try:
                ra_deg, dec_deg = _resolve_target(args.target)
                slug = _make_output_slug(args.target, ra_deg, dec_deg)
            except ValueError as e:
                log.error(str(e))
                return 1
        else:
            if not args.ra_deg or not args.dec_deg:
                log.error("Either --target or both --ra-deg and --dec-deg are required")
                return 1
            ra_deg, dec_deg = args.ra_deg, args.dec_deg
            slug = _make_output_slug(None, ra_deg, dec_deg)
        
        # Resolve width/height from --size-deg or explicit overrides
        width_deg = args.width_deg if args.width_deg is not None else args.size_deg
        height_deg = args.height_deg if args.height_deg is not None else args.size_deg
        
        # Auto-generate output paths
        output_fits = args.output_fits or Path(f"data/output/{slug}_{width_deg}deg.fits")
        manifest_json = args.manifest_json or Path(f"data/output/{slug}_{width_deg}deg_manifest.json")
        
        output_fits.parent.mkdir(parents=True, exist_ok=True)
        manifest_json.parent.mkdir(parents=True, exist_ok=True)
        
        config = BuildConfig(
            region=Region(
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                width_deg=width_deg,
                height_deg=height_deg,
            ),
            as_of_jd=parse_cli_date_jd(args.as_of_date),
            earliest_jd=parse_cli_date_jd(args.earliest_date),
            session_root=args.session_root,
            output_fits=output_fits,
            epoch_fits=args.epoch_fits,
            manifest_json=manifest_json,
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
            ra_deg,
            dec_deg,
            width_deg,
            height_deg,
            args.binning,
        )
        manifest = build_mosaic(config)
        log.info(
            "Build complete: %d plate(s) -> %s",
            len(manifest.get("plates", [])),
            output_fits,
        )
    elif args.command == "plate-photos":
        # Resolve target name if provided
        if args.target and (args.ra_deg or args.dec_deg):
            log.error("--target and --ra-deg/--dec-deg are mutually exclusive")
            return 1
        
        if args.target:
            try:
                ra_deg, dec_deg = _resolve_target(args.target)
                slug = _make_output_slug(args.target, ra_deg, dec_deg)
            except ValueError as e:
                log.error(str(e))
                return 1
        else:
            if not args.ra_deg or not args.dec_deg:
                log.error("Either --target or both --ra-deg and --dec-deg are required")
                return 1
            ra_deg, dec_deg = args.ra_deg, args.dec_deg
            slug = _make_output_slug(None, ra_deg, dec_deg)
        
        # Resolve width/height from --size-deg or explicit overrides
        width_deg = args.width_deg if args.width_deg is not None else args.size_deg
        height_deg = args.height_deg if args.height_deg is not None else args.size_deg
        
        # Auto-generate manifest path
        manifest_json = args.manifest_json or Path(f"data/output/{slug}_{width_deg}deg_photos_manifest.json")
        manifest_json.parent.mkdir(parents=True, exist_ok=True)
        
        config = PlatePhotoConfig(
            region=Region(
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                width_deg=width_deg,
                height_deg=height_deg,
            ),
            as_of_jd=parse_cli_date_jd(args.as_of_date),
            earliest_jd=parse_cli_date_jd(args.earliest_date),
            query_step_deg=args.query_step_deg,
            api_base=args.api_base,
            api_key=os.getenv(args.api_key_env),
            allow_multi_solution_plates=args.allow_multi_solution_plates,
            max_plates=args.max_plates,
            output_dir=args.output_dir,
            manifest_json=manifest_json,
            overwrite=args.overwrite,
        )
        log.info(
            "Starting DASCH plate-photo download: center=(%.6f, %.6f) size=%.3fx%.3f deg",
            ra_deg,
            dec_deg,
            width_deg,
            height_deg,
        )
        manifest = discover_and_download_plate_photos(config)
        log.info(
            "Photo download complete: %d/%d plate(s) have local photos -> %s",
            manifest.get("n_downloaded", 0),
            manifest.get("n_candidates", 0),
            manifest_json,
        )
    elif args.command == "wtml":
        # Resolve target name if provided
        if args.target and (args.ra_deg or args.dec_deg):
            log.error("--target and --ra-deg/--dec-deg are mutually exclusive")
            return 1
        
        region = None
        slug = "wtml"
        
        if args.photo_manifest_json is None:
            # Need to resolve coordinates for discovery
            if args.target:
                try:
                    ra_deg, dec_deg = _resolve_target(args.target)
                    slug = _make_output_slug(args.target, ra_deg, dec_deg)
                except ValueError as e:
                    log.error(str(e))
                    return 1
            else:
                if not args.ra_deg or not args.dec_deg:
                    log.error("Either --target or both --ra-deg and --dec-deg are required for discovery")
                    return 1
                ra_deg, dec_deg = args.ra_deg, args.dec_deg
                slug = _make_output_slug(None, ra_deg, dec_deg)
            
            # Resolve width/height from --size-deg or explicit overrides
            width_deg = args.width_deg if args.width_deg is not None else args.size_deg
            height_deg = args.height_deg if args.height_deg is not None else args.size_deg
            
            region = Region(
                ra_deg=ra_deg,
                dec_deg=dec_deg,
                width_deg=width_deg,
                height_deg=height_deg,
            )

        # Auto-generate output paths
        output_wtml = args.output_wtml or Path(f"data/output/{slug}_{args.size_deg}deg.wtml")
        output_json = args.output_json or Path(f"data/output/{slug}_{args.size_deg}deg_report.json")
        
        output_wtml.parent.mkdir(parents=True, exist_ok=True)
        output_json.parent.mkdir(parents=True, exist_ok=True)

        config = WtmlBuildConfig(
            output_wtml=output_wtml,
            output_json=output_json,
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
            s3_region=args.s3_region,
        )
        log.info("Starting WTML pipeline")
        result = build_wtml(config)
        log.info(
            "WTML complete: paired=%d rejected=%d -> %s",
            result.get("n_paired", 0),
            result.get("n_rejected", 0),
            output_wtml,
        )
    elif args.command == "wtml-solve-plate":
        plate_id = args.plate_id.strip().lower()
        output_wtml = args.output_wtml or Path(f"data/output/{plate_id}_solve.wtml")
        output_json = args.output_json or Path(f"data/output/{plate_id}_solve_report.json")
        output_wtml.parent.mkdir(parents=True, exist_ok=True)
        output_json.parent.mkdir(parents=True, exist_ok=True)

        if (args.ra_hint_deg is None) != (args.dec_hint_deg is None):
            log.error("--ra-hint-deg and --dec-hint-deg must be supplied together")
            return 1
        if (args.scale_low_arcsec_per_pix is None) != (args.scale_high_arcsec_per_pix is None):
            log.error("--scale-low-arcsec-per-pix and --scale-high-arcsec-per-pix must be supplied together")
            return 1

        astrometry_api_key = args.astrometry_api_key or os.getenv(args.astrometry_api_key_env)
        if not astrometry_api_key:
            log.error("Astrometry.net API key is required via --astrometry-api-key or --astrometry-api-key-env")
            return 1

        config = WtmlPlateSolveConfig(
            plate_id=plate_id,
            output_wtml=output_wtml,
            output_json=output_json,
            photo_output_dir=args.photo_output_dir,
            astrometry_api_key=astrometry_api_key,
            api_base=args.api_base,
            api_key=os.getenv(args.api_key_env),
            overwrite=args.overwrite,
            s3_region=args.s3_region,
            astrometry_api_base=args.astrometry_api_base,
            solver_work_dir=args.solver_work_dir,
            solver_timeout_sec=args.solver_timeout_sec,
            solver_upload_timeout_sec=args.solver_upload_timeout_sec,
            solver_poll_interval_sec=args.solver_poll_interval_sec,
            ra_hint_deg=args.ra_hint_deg,
            dec_hint_deg=args.dec_hint_deg,
            radius_hint_deg=args.radius_hint_deg,
            scale_low_arcsec_per_pix=args.scale_low_arcsec_per_pix,
            scale_high_arcsec_per_pix=args.scale_high_arcsec_per_pix,
        )
        log.info("Starting WTML plate solve for %s", plate_id)
        result = build_wtml_from_plate_solve(config)
        log.info(
            "WTML plate solve complete: plate=%s placement_center=(%.6f, %.6f) -> %s",
            plate_id,
            float(result.get("skyimage_placement", {}).get("center_x", 0.0)),
            float(result.get("skyimage_placement", {}).get("center_y", 0.0)),
            output_wtml,
        )
    else:  # pragma: no cover
        parser.error(f"unsupported command: {args.command}")
    return 0