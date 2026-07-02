"""Build WorldWide Telescope webclient links for single DASCH plates.

WWT's webclient can load a single image directly via its ShowImage.aspx endpoint:
    https://www.worldwidetelescope.org/webclient/?wtml=<url-encoded ShowImage.aspx URL>

This module builds that link.
"""
from __future__ import annotations

from typing import Any
from urllib.parse import urlencode

from astropy.wcs import WCS
from wwt_data_formats.imageset import ImageSet as WwtImageSet

SHOWIMAGE_URL = "http://www.worldwidetelescope.org/wwtweb/ShowImage.aspx"
WEBCLIENT_URL = "https://www.worldwidetelescope.org/webclient/"

DEFAULT_CREDITS = "DASCH / Harvard College Observatory"
DEFAULT_CREDITS_URL = "https://dasch.cfa.harvard.edu/"


def derive_showimage_placement(
    wcs_header: dict[str, Any],
    width_px: int,
    height_px: int,
) -> dict[str, float | bool]:
    """Derive ShowImage.aspx placement fields (scale, rotation, ra/dec, parity, ref pixel) from a WCS header.

    wwt_data_formats' rotation convention is offset by 180 degrees from what
    ShowImage.aspx expects (found empirically).

    Delegates the CD-matrix decomposition (scale, rotation, parity) to
    wwt_data_formats rather than hand-rolling the trig, since getting the
    mirrored-image (reverseparity) case right requires care.
    """
    imgset = WwtImageSet()
    imgset.tile_levels = 0
    imgset.set_position_from_wcs(wcs_header, width=width_px, height=height_px)

    center_x_px = width_px / 2.0
    center_y_px = height_px / 2.0
    # Plain-array WCS math (no SkyCoord/Quantity) -- avoids astropy's notoriously
    # loose Quantity/SkyCoord type stubs; we just want two degree values.
    ra_deg, dec_deg = WCS(wcs_header).all_pix2world(center_x_px, center_y_px, 0)

    return {
        "scale_arcsec_per_px": float(imgset.base_degrees_per_tile) * 3600.0,
        "rotation_deg": (float(imgset.rotation_deg) + 180.0) % 360.0,
        "ra_deg": float(ra_deg),
        "dec_deg": float(dec_deg),
        "ref_pixel_x": center_x_px,
        "ref_pixel_y": center_y_px,
        "reverseparity": bool(imgset.bottoms_up),
    }


def build_wtml_link(
    name: str,
    imageurl: str,
    ra: float,
    dec: float,
    x: float,
    y: float,
    rotation: float,
    scale: float,
    reverseparity: bool,
    thumb: str = "",
    credits: str = DEFAULT_CREDITS,
    creditsurl: str = DEFAULT_CREDITS_URL,
) -> str:
    """Build a WWT webclient link that opens one plate image at its sky position.

    Args:
        name: Plate identifier, shown as the image name in WWT.
        imageurl: Publicly fetchable URL of the plate JPG.
        ra, dec: Sky position of the reference pixel (x, y), in degrees.
        x, y: Reference pixel coordinates within the full-resolution image
            that (ra, dec) corresponds to -- NOT the image's pixel dimensions.
        rotation: Rotation on sky, in degrees (East from North).
        scale: Pixel scale, in arcsec/pixel.
        reverseparity: Whether the image is mirrored
        thumb: Optional thumbnail URL.
        credits, creditsurl: Attribution shown in WWT.
    """
    show_image_params = {
        "reverseparity": str(reverseparity),
        "scale": f"{scale:.5f}",
        "name": name,
        "imageurl": imageurl,
        "credits": credits,
        "creditsUrl": creditsurl,
        "ra": f"{ra:.9f}",
        "y": str(y),
        "x": str(x),
        "rotation": f"{rotation:.4f}",
        "dec": f"{dec:.9f}",
        "thumb": thumb,
        "wtml": "true",
    }
    show_image_url = f"{SHOWIMAGE_URL}?{urlencode(show_image_params)}"
    return f"{WEBCLIENT_URL}?{urlencode({'wtml': show_image_url})}"
