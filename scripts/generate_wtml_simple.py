from __future__ import annotations

import json
from pathlib import Path
import xml.etree.ElementTree as ET

from dasch_sky_mosaic.wtml import (
    _load_plate_wcs,
    _prepare_photo_for_wcs,
    _derive_skyimage_from_wcs,
    _make_imageset_xml,
)


def find_file(root: Path, prefix: str) -> Path | None:
    for p in root.iterdir():
        if p.name.lower().startswith(prefix.lower()):
            return p
    return None


def main(plate_id: str) -> int:
    plate_id = plate_id.lower()
    photo_dir = Path("data/cache/dasch_session/plate_photos")
    mosaic_dir = Path("data/cache/dasch_session/mosaic_package")
    photo = find_file(photo_dir, plate_id)
    mosaic = find_file(mosaic_dir, plate_id)
    if photo is None or mosaic is None:
        print("Missing photo or mosaic for", plate_id)
        return 1

    header, wcs, shape, scale_deg, fits_data = _load_plate_wcs(mosaic)
    # Force fallback placement by not providing fits_data to avoid astroalign
    image_path, image_w, image_h, placement_header, alignment_meta = _prepare_photo_for_wcs(
        photo.resolve(), wcs, shape, fits_data=None
    )
    placement = _derive_skyimage_from_wcs(placement_header, image_w, image_h)

    root = ET.Element(
        "Folder",
        {
            "Name": "DASCH Photo WTML",
            "Group": "Explorer",
            "Type": "Sky",
            "Searchable": "True",
        },
    )

    place = _make_imageset_xml(plate_id, image_path.resolve().as_uri(), placement, image_w, image_h)
    root.append(place)

    out_wtml = Path("data/output") / f"{plate_id}.wtml"
    out_json = Path("data/output") / f"{plate_id}_wtml.json"
    out_wtml.parent.mkdir(parents=True, exist_ok=True)
    xml_text = ET.tostring(root, encoding="unicode")
    out_wtml.write_text("<?xml version='1.0' encoding='UTF-8'?>\n" + xml_text + "\n", encoding="utf-8")

    result = {
        "workflow": "wtml-build-simple",
        "n_paired": 1,
        "paired": [
            {
                "plate_id": plate_id,
                "local_photo_path": str(photo),
                "aligned_photo_path": str(image_path),
                "wcs_mosaic_path": str(mosaic),
                "alignment": alignment_meta,
                "skyimage_placement": placement,
            }
        ],
        "output_wtml": str(out_wtml),
    }
    out_json.write_text(json.dumps(result, indent=2), encoding="utf-8")
    print(out_wtml, out_json)
    return 0


if __name__ == "__main__":
    import sys

    if len(sys.argv) < 2:
        print("Usage: generate_wtml_simple.py <plate_id>")
        raise SystemExit(2)
    raise SystemExit(main(sys.argv[1]))
