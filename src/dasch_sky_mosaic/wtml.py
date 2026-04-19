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
from wwt_data_formats.imageset import ImageSet as WwtImageSet

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


def _load_plate_wcs(mosaic_path: Path) -> tuple[Any, WCS, tuple[int, int], float, "Any"]:
    with fits.open(mosaic_path, memmap=True) as hdul:
        hdu = hdul[0]  # type: ignore[index]
        header = hdu.header.copy()  # type: ignore[attr-defined]
        data = hdu.data  # type: ignore[attr-defined]
        if data is None:
            raise RuntimeError(f"missing image data in {mosaic_path}")
        shape = (int(data.shape[0]), int(data.shape[1]))
        # Make a copy of data to keep it in memory after file is closed
        data_copy = data.copy()
    wcs = WCS(header)
    wcs_header = wcs.to_header(relax=True)
    scale_deg = float(abs(proj_plane_pixel_scales(wcs).mean()))
    return wcs_header, wcs, shape, scale_deg, data_copy


def _next_pow2(x: int) -> int:
    return 1 if x <= 1 else 1 << (x - 1).bit_length()


def _astroalign_find_transform(
    alignment_source: "Any",
    fits_data: "Any",
) -> tuple["Any" | None, int, str]:
    """Find a similarity transform from a plate-photo image to FITS pixels.

    Returns (transform, n_matches, mode_label).
    """
    try:
        import astroalign as aa
        import numpy as np
    except ImportError:
        LOG.warning("Astroalign alignment unavailable; missing dependencies.")
        return None, 0, "unavailable"

    def _normalize(img: "Any") -> "Any":
        p1, p99 = np.percentile(img, (1.0, 99.0))
        if p99 <= p1:
            return np.zeros_like(img, dtype=np.float32)
        img = np.clip(img, p1, p99)
        return ((img - p1) / (p99 - p1)).astype(np.float32)

    try:
        src_rgb = np.array(alignment_source, dtype=np.float32)
        src_gray = _normalize(src_rgb[:, :, 0])
        fits_norm = _normalize(np.array(fits_data, dtype=np.float32))

        best_transform: Any = None
        best_label = ""
        best_matches = 0

        # Inversion is tested on the plate-photo source, not on FITS.
        for label, source in (("native", src_gray), ("inverted", 1.0 - src_gray)):
            try:
                transform, (src_pts, dst_pts) = aa.find_transform(
                    source,
                    fits_norm,
                    max_control_points=150,
                    detection_sigma=4,
                    min_area=3,
                )
            except Exception:
                continue

            n_matches = int(len(src_pts))
            if n_matches > best_matches:
                best_matches = n_matches
                best_transform = transform
                best_label = label

        if best_transform is None or best_matches < 6:
            return None, 0, "none"

        LOG.info(
            "Astroalign alignment: mode=%s matches=%d scale=%.5f rot_deg=%.3f tx=%.2f ty=%.2f",
            best_label,
            best_matches,
            float(best_transform.scale),
            float(np.degrees(best_transform.rotation)),
            float(best_transform.translation[0]),
            float(best_transform.translation[1]),
        )
        return best_transform, best_matches, best_label
    except Exception:
        LOG.exception("Astroalign alignment failed; using unaligned crop.")
        return None, 0, "error"


def _compose_wcs_header_from_affine(reference_wcs: WCS, affine_a: "Any", affine_t: "Any") -> dict[str, Any]:
    """Compose an image->reference affine transform into a FITS WCS header.

    The affine maps source image x/y pixel coordinates (0-based) into reference
    FITS x/y pixel coordinates (0-based).
    """
    import numpy as np

    if reference_wcs.wcs.has_cd():
        ref_cd = np.array(reference_wcs.wcs.cd, dtype=np.float64)
    else:
        ref_cd = np.array(reference_wcs.wcs.get_pc(), dtype=np.float64) @ np.diag(
            np.array(reference_wcs.wcs.cdelt, dtype=np.float64)
        )

    ones = np.ones(2, dtype=np.float64)
    affine_a = np.array(affine_a, dtype=np.float64)
    affine_t = np.array(affine_t, dtype=np.float64)
    cd_src = ref_cd @ affine_a
    b = affine_t + ones - affine_a @ ones
    crpix_src = np.linalg.solve(affine_a, np.array(reference_wcs.wcs.crpix, dtype=np.float64) - b)

    header = reference_wcs.to_header(relax=True)
    for key in list(header.keys()):
        if key.startswith("PC") or key.startswith("CD"):
            del header[key]
    for key in ("CDELT1", "CDELT2"):
        if key in header:
            del header[key]

    header["CRPIX1"] = float(crpix_src[0])
    header["CRPIX2"] = float(crpix_src[1])
    header["CD1_1"] = float(cd_src[0, 0])
    header["CD1_2"] = float(cd_src[0, 1])
    header["CD2_1"] = float(cd_src[1, 0])
    header["CD2_2"] = float(cd_src[1, 1])
    return dict(header)


def _prepare_photo_for_wcs(
    source_photo: Path,
    reference_wcs: WCS,
    shape: tuple[int, int],
    fits_data: "Any" = None,
) -> tuple[Path, int, int, dict[str, Any], dict[str, Any]]:
    """Map a raw plate photo onto the FITS WCS pixel geometry.

    When FITS data is available, we compare two parity branches only: the source
    image and its horizontal mirror. Astroalign then solves rotation, translation,
    and uniform scale within each branch, and we keep the stronger match.
    """
    try:
        from PIL import Image
        import numpy as np
    except ImportError as exc:  # pragma: no cover
        raise RuntimeError("Pillow is required for WTML tile generation. Install pillow>=10.") from exc

    import numpy as np

    target_h, target_w = shape
    with Image.open(source_photo) as opened:
        src = opened.convert("RGB")
    src_w, src_h = src.size

    def _cover_resize(candidate: "Any") -> tuple["Any", float]:
        c_w, c_h = candidate.size
        scale = max(target_w / max(1, c_w), target_h / max(1, c_h))
        resized_w = max(1, int(round(c_w * scale)))
        resized_h = max(1, int(round(c_h * scale)))
        resized = candidate.resize((resized_w, resized_h), resample=Image.Resampling.BICUBIC)
        return resized, float(scale)

    def _branch_affine(branch: str, width_px: int, scale: float) -> tuple["Any", "Any"]:
        if branch == "source_mirror":
            base_a = np.array([[-1.0, 0.0], [0.0, 1.0]], dtype=np.float64)
            base_t = np.array([width_px - 1.0, 0.0], dtype=np.float64)
        else:
            base_a = np.eye(2, dtype=np.float64)
            base_t = np.zeros(2, dtype=np.float64)
        return scale * base_a, scale * base_t

    if fits_data is not None:
        parity_candidates = [
            ("source", src),
            ("source_mirror", src.transpose(Image.Transpose.FLIP_LEFT_RIGHT)),
        ]
        best_header = None
        best_label = ""
        best_matches = -1
        best_meta: dict[str, Any] | None = None

        for label, candidate in parity_candidates:
            resized, resize_scale = _cover_resize(candidate)
            transform, n_matches, mode_label = _astroalign_find_transform(resized, fits_data)
            if transform is None:
                continue

            branch_a, branch_t = _branch_affine(label, src_w, resize_scale)
            total_a = np.array(transform.params[:2, :2], dtype=np.float64) @ branch_a
            total_t = np.array(transform.params[:2, 2], dtype=np.float64) + np.array(transform.params[:2, :2], dtype=np.float64) @ branch_t
            header = _compose_wcs_header_from_affine(reference_wcs, total_a, total_t)

            if n_matches > best_matches:
                best_header = header
                best_label = f"{label}/{mode_label}"
                best_matches = n_matches
                best_meta = {
                    "parity": label,
                    "mode": mode_label,
                    "matches": n_matches,
                    "transform_scale": float(transform.scale),
                    "transform_rotation_deg": float(np.degrees(transform.rotation)),
                    "transform_translation": [
                        float(transform.translation[0]),
                        float(transform.translation[1]),
                    ],
                }

        if best_header is not None and best_meta is not None:
            LOG.info("Parity pick: %s (matches=%d)", best_label, best_matches)
            return source_photo, src_w, src_h, best_header, best_meta

    resized, resize_scale = _cover_resize(src)
    x0 = max(0, min((resized.width - target_w) // 2, resized.width - target_w))
    y0 = max(0, min((resized.height - target_h) // 2, resized.height - target_h))
    fallback_a = resize_scale * np.eye(2, dtype=np.float64)
    fallback_t = np.array([-x0, -y0], dtype=np.float64)
    fallback_header = _compose_wcs_header_from_affine(reference_wcs, fallback_a, fallback_t)
    fallback_meta = {
        "parity": "source",
        "mode": "fallback",
        "matches": 0,
        "transform_scale": 1.0,
        "transform_rotation_deg": 0.0,
        "transform_translation": [float(-x0), float(-y0)],
    }
    LOG.info("Parity pick fallback: source")
    return source_photo, src_w, src_h, fallback_header, fallback_meta


def _write_aligned_photo(aligned_rgb: "Any", out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    aligned_rgb.save(out_path, format="JPEG", quality=92)
    return out_path


def _derive_skyimage_from_wcs(
    header: dict[str, Any],
    width_px: int,
    height_px: int,
) -> dict[str, float | bool]:
    """Derive WWT SkyImage placement fields from FITS WCS and image dimensions."""
    imgset = WwtImageSet()
    imgset.tile_levels = 0
    imgset.set_position_from_wcs(header, width=width_px, height=height_px)

    return {
        "base_degrees_per_tile": float(imgset.base_degrees_per_tile),
        "bottoms_up": bool(imgset.bottoms_up),
        "center_x": float(imgset.center_x),
        "center_y": float(imgset.center_y),
        "offset_x": float(imgset.offset_x),
        "offset_y": float(imgset.offset_y),
        "rotation_deg": float(imgset.rotation_deg),
        "width_factor": float(imgset.width_factor),
    }


def _make_imageset_xml(
    plate_id: str,
    image_url: str,
    placement: dict[str, float | bool],
    image_width_px: int,
    image_height_px: int,
) -> ET.Element:
    center_ra_deg = float(placement["center_x"])
    center_dec_deg = float(placement["center_y"])
    base_deg_per_px = max(1e-12, float(placement["base_degrees_per_tile"]))
    width_deg = image_width_px * base_deg_per_px
    height_deg = image_height_px * base_deg_per_px

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
            "BaseDegreesPerTile": f"{base_deg_per_px:.12f}",
            "FileType": ".jpg",
            "Projection": "SkyImage",
            "BottomsUp": "True" if bool(placement["bottoms_up"]) else "False",
            "CenterX": f"{center_ra_deg:.8f}",
            "CenterY": f"{center_dec_deg:.8f}",
            "OffsetX": f"{float(placement['offset_x']):.8f}",
            "OffsetY": f"{float(placement['offset_y']):.8f}",
            "Rotation": f"{float(placement['rotation_deg']):.8f}",
            "Sparse": "False",
            "WidthFactor": f"{float(placement['width_factor']):.8f}",
            "QuadTreeMap": "",
        },
    )
    ET.SubElement(imageset, "ThumbnailUrl")
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

        header, wcs, shape, scale_deg, fits_data = _load_plate_wcs(mosaic_path)
        image_path, image_w, image_h, placement_header, alignment_meta = _prepare_photo_for_wcs(
            photo_file.resolve(),
            wcs,
            shape,
            fits_data=fits_data,
        )
        placement = _derive_skyimage_from_wcs(placement_header, image_w, image_h)
        image_uri = image_path.resolve().as_uri()
        center_ra_deg = float(wcs.wcs.crval[0])
        center_dec_deg = float(wcs.wcs.crval[1])

        root.append(
            _make_imageset_xml(
                plate_id=plate_id,
                image_url=image_uri,
                placement=placement,
                image_width_px=image_w,
                image_height_px=image_h,
            )
        )

        paired_rows.append(
            {
                "plate_id": plate_id,
                "local_photo_path": photo_path,
                "aligned_photo_path": str(image_path),
                "wcs_mosaic_path": str(mosaic_path),
                "center_ra_deg": center_ra_deg,
                "center_dec_deg": center_dec_deg,
                "pixel_scale_deg": scale_deg,
                "shape": [shape[0], shape[1]],
                "alignment": alignment_meta,
                "skyimage_placement": placement,
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
