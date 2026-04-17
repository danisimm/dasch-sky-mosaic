from __future__ import annotations

import json
import logging
import math
import xml.etree.ElementTree as ET
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales

from dasch_sky_mosaic.fetch import Region, download_mosaic_paths, parse_cli_date_jd
from dasch_sky_mosaic.plate_photos import PlatePhotoConfig, discover_and_download_plate_photos

LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class WtmlBuildConfig:
    output_wtml: Path
    output_json: Path
    photo_manifest_json: Path | None
    region: Region | None
    as_of_date: str | None
    earliest_date: str | None
    query_step_deg: float
    api_base: str
    api_key: str | None
    allow_multi_solution_plates: bool
    max_plates: int | None
    session_root: Path
    photo_output_dir: Path
    overwrite: bool


def _require_region(config: WtmlBuildConfig) -> Region:
    if config.region is None:
        raise ValueError("region is required when photo_manifest_json is not provided")
    return config.region


def _discover_or_load_photos(config: WtmlBuildConfig) -> dict[str, Any]:
    if config.photo_manifest_json is not None:
        return json.loads(config.photo_manifest_json.read_text(encoding="utf-8"))

    region = _require_region(config)
    auto_manifest_path = config.output_json.with_name(config.output_json.stem + "_photos.json")
    photo_cfg = PlatePhotoConfig(
        region=region,
        as_of_jd=parse_cli_date_jd(config.as_of_date),
        earliest_jd=parse_cli_date_jd(config.earliest_date),
        query_step_deg=config.query_step_deg,
        api_base=config.api_base,
        api_key=config.api_key,
        allow_multi_solution_plates=config.allow_multi_solution_plates,
        max_plates=config.max_plates,
        output_dir=config.photo_output_dir,
        manifest_json=auto_manifest_path,
        overwrite=config.overwrite,
    )
    LOG.info("No photo manifest supplied; running fresh discovery and photo download")
    return discover_and_download_plate_photos(photo_cfg)


def _load_plate_wcs(mosaic_path: Path) -> tuple[WCS, tuple[int, int], float]:
    with fits.open(mosaic_path, memmap=True) as hdul:
        hdu = hdul[0]  # type: ignore[index]
        header = hdu.header.copy()  # type: ignore[attr-defined]
        data = hdu.data  # type: ignore[attr-defined]
        if data is None:
            raise RuntimeError(f"missing image data in {mosaic_path}")
        shape = (int(data.shape[0]), int(data.shape[1]))
    wcs = WCS(header)
    scale_deg = float(abs(proj_plane_pixel_scales(wcs).mean()))
    return wcs, shape, scale_deg


def _next_pow2(x: int) -> int:
    return 1 if x <= 1 else 1 << (x - 1).bit_length()


def _prepare_photo_for_wcs(source_photo: Path, shape: tuple[int, int]) -> "Any":
    """Map a raw plate photo onto the FITS WCS pixel geometry.

    We do a center-preserving cover fit, with optional 90-degree rotation if that
    better matches the target aspect ratio. This keeps a direct photo workflow
    while ensuring tile pixels correspond to the same sky footprint as the WCS.
    """
    try:
        from PIL import Image
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Pillow is required for WTML tile generation. Install pillow>=10.") from exc

    target_h, target_w = shape
    with Image.open(source_photo) as opened:
        src = opened.convert("RGB")

    def _score(candidate: "Any") -> float:
        c_w, c_h = candidate.size
        src_aspect = c_w / max(1.0, c_h)
        dst_aspect = target_w / max(1.0, target_h)
        return abs(math.log(max(1e-12, src_aspect / dst_aspect)))

    candidates = [src, src.transpose(Image.Transpose.ROTATE_90)]
    best = min(candidates, key=_score)
    b_w, b_h = best.size
    scale = max(target_w / max(1, b_w), target_h / max(1, b_h))
    resized_w = max(1, int(round(b_w * scale)))
    resized_h = max(1, int(round(b_h * scale)))
    resized = best.resize((resized_w, resized_h), resample=Image.Resampling.BICUBIC)

    x0 = max(0, (resized_w - target_w) // 2)
    y0 = max(0, (resized_h - target_h) // 2)
    return resized.crop((x0, y0, x0 + target_w, y0 + target_h))


def _write_aligned_photo(aligned_rgb: "Any", out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    aligned_rgb.save(out_path, format="JPEG", quality=92)
    return out_path


def _make_imageset_xml(
    plate_id: str,
    image_url: str,
    center_ra_deg: float,
    center_dec_deg: float,
    shape: tuple[int, int],
    scale_deg: float,
) -> ET.Element:
    naxis2, naxis1 = shape
    width_deg = max(1e-9, naxis1 * scale_deg)
    height_deg = max(1e-9, naxis2 * scale_deg)
    width_factor = naxis1 / max(1.0, naxis2)
    # SkyImage anchors from an image origin unless an offset is supplied; shift by
    # half-extents so the image center lands on (CenterX, CenterY).
    offset_x_deg = width_deg * 0.5
    offset_y_deg = height_deg * 0.5

    place = ET.Element(
        "Place",
        {
            "Name": plate_id,
            "DataSetType": "Sky",
            "RA": f"{center_ra_deg / 15.0:.8f}",
            "Dec": f"{center_dec_deg:.8f}",
            "ZoomLevel": f"{max(width_deg, height_deg) * 6.0:.6f}",
        },
    )
    fg = ET.SubElement(place, "ForegroundImageSet")
    imageset = ET.SubElement(
        fg,
        "ImageSet",
        {
            "Generic": "False",
            "DataSetType": "Sky",
            "BandPass": "Visible",
            "Name": f"DASCH {plate_id}",
            "Url": image_url,
            "TileLevels": "0",
            "BaseTileLevel": "0",
            "BaseDegreesPerTile": f"{height_deg:.12f}",
            "FileType": ".jpg",
            "Projection": "SkyImage",
            "BottomsUp": "False",
            "CenterX": f"{center_ra_deg:.8f}",
            "CenterY": f"{center_dec_deg:.8f}",
            "OffsetX": f"{offset_x_deg:.8f}",
            "OffsetY": f"{offset_y_deg:.8f}",
            "Rotation": "0",
            "Sparse": "True",
            "WidthFactor": f"{width_factor:.8f}",
        },
    )
    ET.SubElement(imageset, "Credits").text = "DASCH / Harvard College Observatory"
    ET.SubElement(imageset, "CreditsUrl").text = "https://dasch.cfa.harvard.edu/"
    ET.SubElement(imageset, "Description").text = f"Single-plate WTML entry for {plate_id}"
    return place


def build_wtml(config: WtmlBuildConfig) -> dict[str, Any]:
    photos_manifest = _discover_or_load_photos(config)
    rows = [r for r in photos_manifest.get("photos", []) if r.get("status") in {"downloaded", "cached"}]
    if not rows:
        raise RuntimeError("no successful photo rows available")

    plate_ids = [str(r["plate_id"]).lower() for r in rows]
    mosaic_paths = download_mosaic_paths(
        session_root=config.session_root,
        plate_ids=plate_ids,
        binning=16,
        api_base=config.api_base,
        api_key=config.api_key,
    )

    root = ET.Element(
        "Folder",
        {
            "Name": "DASCH Photo WTML",
            "Group": "Explorer",
            "Type": "Sky",
            "Searchable": "True",
        },
    )

    paired_rows: list[dict[str, Any]] = []
    rejected_rows: list[dict[str, Any]] = []

    for row in rows:
        plate_id = str(row["plate_id"]).lower()
        photo_path = str(row.get("local_photo_path", ""))
        photo_file = Path(photo_path)
        if not photo_file.exists():
            rejected_rows.append(
                {
                    "plate_id": plate_id,
                    "reason": "missing_local_photo_file",
                    "local_photo_path": photo_path,
                }
            )
            continue

        mosaic_path = mosaic_paths.get(plate_id)
        if mosaic_path is None:
            rejected_rows.append(
                {
                    "plate_id": plate_id,
                    "reason": "missing_mosaic_wcs",
                    "local_photo_path": photo_path,
                }
            )
            continue

        wcs, shape, scale_deg = _load_plate_wcs(mosaic_path)
        aligned_photo = _prepare_photo_for_wcs(photo_file.resolve(), shape)
        aligned_dir = config.output_wtml.parent / f"{config.output_wtml.stem}_images"
        aligned_path = _write_aligned_photo(aligned_photo, aligned_dir / f"{plate_id}.jpg")
        image_uri = aligned_path.resolve().as_uri()
        center_ra_deg = float(wcs.wcs.crval[0])
        center_dec_deg = float(wcs.wcs.crval[1])

        root.append(
            _make_imageset_xml(
                plate_id=plate_id,
                image_url=image_uri,
                center_ra_deg=center_ra_deg,
                center_dec_deg=center_dec_deg,
                shape=shape,
                scale_deg=scale_deg,
            )
        )

        paired_rows.append(
            {
                "plate_id": plate_id,
                "local_photo_path": photo_path,
                "aligned_photo_path": str(aligned_path),
                "wcs_mosaic_path": str(mosaic_path),
                "center_ra_deg": center_ra_deg,
                "center_dec_deg": center_dec_deg,
                "pixel_scale_deg": scale_deg,
                "shape": [shape[0], shape[1]],
            }
        )

    config.output_wtml.parent.mkdir(parents=True, exist_ok=True)
    config.output_json.parent.mkdir(parents=True, exist_ok=True)
    if (config.output_wtml.exists() or config.output_json.exists()) and not config.overwrite:
        raise FileExistsError("refusing to overwrite output artifacts; rerun with --overwrite")

    xml_text = ET.tostring(root, encoding="unicode")
    config.output_wtml.write_text("<?xml version='1.0' encoding='UTF-8'?>\n" + xml_text + "\n", encoding="utf-8")

    result: dict[str, Any] = {
        "workflow": "wtml-build",
        "source_photo_manifest": str(config.photo_manifest_json) if config.photo_manifest_json else "generated-in-run",
        "n_photo_rows": len(rows),
        "n_paired": len(paired_rows),
        "n_rejected": len(rejected_rows),
        "paired": paired_rows,
        "rejected": rejected_rows,
        "output_wtml": str(config.output_wtml),
    }
    config.output_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result
