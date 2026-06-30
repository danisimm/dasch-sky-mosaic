"""Pure WTML XML generation.

Takes pre-computed WCS placement data and produces WTML output compatible
with WorldWide Telescope. No downloads or solving happen here.
"""
from __future__ import annotations

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from wwt_data_formats.imageset import ImageSet as WwtImageSet


def derive_skyimage_placement(
    wcs_header: dict[str, Any],
    width_px: int,
    height_px: int,
) -> dict[str, float | bool]:
    """Derive WWT SkyImage placement fields from a FITS WCS header and image dimensions."""
    imgset = WwtImageSet()
    imgset.tile_levels = 0
    imgset.set_position_from_wcs(wcs_header, width=width_px, height=height_px)
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


def make_place_element(
    plate_id: str,
    image_url: str,
    placement: dict[str, float | bool],
    image_width_px: int,
    image_height_px: int,
) -> ET.Element:
    """Build a WWT <Place> XML element containing an <ImageSet> for one plate."""
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
            "Opacity": "100",
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


def write_wtml(
    plates: list[dict[str, Any]],
    output_path: Path,
    folder_name: str = "DASCH Plates",
) -> None:
    """Write a WTML file from a list of solved plate records.

    Each record in `plates` must have:
      - plate_id: str
      - image_url: str
      - wcs_header: dict
      - image_width_px: int
      - image_height_px: int
    """
    root = ET.Element(
        "Folder",
        {
            "Name": folder_name,
            "Group": "Explorer",
            "Type": "Sky",
            "Searchable": "True",
        },
    )
    for p in plates:
        placement = derive_skyimage_placement(
            p["wcs_header"], int(p["image_width_px"]), int(p["image_height_px"])
        )
        root.append(make_place_element(
            plate_id=str(p["plate_id"]),
            image_url=str(p["image_url"]),
            placement=placement,
            image_width_px=int(p["image_width_px"]),
            image_height_px=int(p["image_height_px"]),
        ))

    output_path.parent.mkdir(parents=True, exist_ok=True)
    xml_text = ET.tostring(root, encoding="unicode")
    output_path.write_text(
        "<?xml version='1.0' encoding='UTF-8'?>\n" + xml_text + "\n",
        encoding="utf-8",
    )
