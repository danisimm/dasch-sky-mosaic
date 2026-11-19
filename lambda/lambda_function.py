"""AWS Lambda: open a DASCH plate JPG in the WorldWide Telescope webclient.

Requests:
  <function-url>/?plate=b02312              -> 302 redirect to the WWT webclient (the website link)
  <function-url>/?plate=b02312&format=wtml  -> WTML for that plate (what the webclient loads)
  <function-url>/?plate=b02312&format=json  -> JSON with the WTML and link (for testing)

Steps:
  1. If the plate was solved before, load the saved WCS from S3.
  2. Otherwise read the plate JPG and its FITS mosaic from S3 into memory, align
     them with astroalign to get the JPG's WCS, and save it to S3.
  3. Redirect to WWT with ?wtml=<this function>?plate=...&format=wtml. The WTML is
     generated on each request, with a fresh presigned URL for the JPG, so the public
     link contains no signed URL and never goes stale.

S3 access comes from the Lambda's IAM role; there are no credentials in this file.
"""
from __future__ import annotations

import io
import json
import logging
import multiprocessing as mp
import os
import re
import warnings
from typing import Any
from urllib.parse import urlencode
from xml.etree import ElementTree as ET

import boto3
import numpy as np
from botocore.exceptions import ClientError
from astropy.io import fits
from astropy.wcs import FITSFixedWarning, WCS
from PIL import Image, ImageOps
from wwt_data_formats.imageset import ImageSet as WwtImageSet

LOG = logging.getLogger()
LOG.setLevel(logging.INFO)

PLATES_BUCKET = os.environ.get("PLATES_BUCKET", "dasch-prod-user")
RESULTS_BUCKET = os.environ.get("RESULTS_BUCKET", PLATES_BUCKET)
RESULTS_PREFIX = os.environ.get("RESULTS_PREFIX", "wwt-wcs/")
IMAGE_URL_EXPIRY_S = int(os.environ.get("IMAGE_URL_EXPIRY_S", "43200"))  # 12 h

PLATE_ID_RE = re.compile(r"^[a-z]{1,5}\d{5}$")
ASTROALIGN_TIMEOUT_S = 5.0

WEBCLIENT_URL = "https://www.worldwidetelescope.org/webclient/"
CREDITS = "DASCH / Harvard College Observatory"
CREDITS_URL = "https://dasch.cfa.harvard.edu/"

s3 = boto3.client("s3", region_name="us-east-1")


# --- Lambda entry point -------------------------------------------------------

def lambda_handler(event: dict[str, Any], context: Any) -> dict[str, Any]:
    params = event.get("queryStringParameters") or {}
    plate_id = str(params.get("plate", "")).strip().lower()
    output = params.get("format", "redirect")  # redirect | wtml | json

    if not PLATE_ID_RE.match(plate_id):
        return _text(400, "Missing or invalid ?plate= (expected something like b02312).")

    try:
        result = load_saved_result(plate_id)
        if result is None:
            LOG.info("Solving %s", plate_id)
            result = solve_plate(plate_id)
            save_result(plate_id, result)
        else:
            LOG.info("Using saved solution for %s", plate_id)
    except FileNotFoundError as exc:
        return _text(404, str(exc))
    except RuntimeError as exc:
        LOG.warning("No solution for %s: %s", plate_id, exc)
        return _text(422, f"Could not place plate {plate_id} on the sky: {exc}")

    wwt_url = build_wwt_link(event, plate_id)

    if output == "wtml":
        # What the WWT webclient fetches. CORS header so worldwidetelescope.org may read it.
        return {
            "statusCode": 200,
            "headers": {"Content-Type": "text/xml", "Access-Control-Allow-Origin": "*"},
            "body": build_wtml(plate_id, result),
        }
    if output == "json":
        return {
            "statusCode": 200,
            "headers": {"Content-Type": "application/json"},
            "body": json.dumps(
                {"plate": plate_id, "wwt_url": wwt_url, "wtml": build_wtml(plate_id, result), "result": result}
            ),
        }
    if wwt_url is None:
        return _text(500, "Can't work out this function's own URL to build the WWT link.")
    return {"statusCode": 302, "headers": {"Location": wwt_url}, "body": ""}


def _text(status: int, message: str) -> dict[str, Any]:
    return {"statusCode": status, "headers": {"Content-Type": "text/plain"}, "body": message}


# --- S3 ----------------------------------------------------------------------

def _result_key(plate_id: str) -> str:
    return f"{RESULTS_PREFIX}{plate_id}.json"


def _photo_key(plate_id: str) -> str:
    return f"plates/{plate_id}/{plate_id}_pphoto_all.jpg"


def load_saved_result(plate_id: str) -> dict[str, Any] | None:
    """Return the saved solution, or None if there isn't one (or it can't be read)."""
    try:
        obj = s3.get_object(Bucket=RESULTS_BUCKET, Key=_result_key(plate_id))
    except s3.exceptions.NoSuchKey:
        return None
    except ClientError as exc:
        # Usually a permissions problem; solving still works, so don't fail the click.
        LOG.warning("Could not read saved result for %s: %s", plate_id, exc)
        return None
    return json.loads(obj["Body"].read())


def save_result(plate_id: str, result: dict[str, Any]) -> None:
    """Save the solution so later clicks skip solving. A failure only logs a warning."""
    try:
        s3.put_object(
            Bucket=RESULTS_BUCKET,
            Key=_result_key(plate_id),
            Body=json.dumps(result, indent=2).encode(),
            ContentType="application/json",
        )
    except ClientError as exc:
        LOG.warning("Could not save result for %s: %s", plate_id, exc)


def _read_s3(key: str) -> bytes:
    try:
        return s3.get_object(Bucket=PLATES_BUCKET, Key=key)["Body"].read()
    except s3.exceptions.NoSuchKey:
        raise FileNotFoundError(f"Not found in S3: {key}") from None


def read_plate_photo(plate_id: str) -> Image.Image:
    data = _read_s3(_photo_key(plate_id))
    with Image.open(io.BytesIO(data)) as opened:
        # Browsers/WWT rotate a displayed image by its EXIF orientation tag; PIL
        # doesn't. Normalize so we solve the orientation that will be shown.
        return ImageOps.exif_transpose(opened).convert("L")


def read_plate_fits(plate_id: str) -> tuple[WCS, np.ndarray]:
    """Return (wcs, image) from the plate's binning-16 FITS mosaic (.fit.fz)."""
    prefix = f"plates/{plate_id}/{plate_id}_mosaic_"
    listing = s3.list_objects_v2(Bucket=PLATES_BUCKET, Prefix=prefix)
    keys = [o["Key"] for o in listing.get("Contents", []) if o["Key"].endswith("ww.fit.fz")]
    if not keys:
        raise FileNotFoundError(f"No binning-16 FITS mosaic for {plate_id}")

    data = _read_s3(sorted(keys)[0])
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=FITSFixedWarning)
        with fits.open(io.BytesIO(data)) as hdul:
            # .fz files keep the image in a compressed extension, not HDU 0.
            hdu = next(h for h in hdul if h.data is not None and h.data.ndim == 2)
            return WCS(hdu.header), np.array(hdu.data, dtype=np.float32)


# --- Alignment (from src/dasch_sky_mosaic/wcs_solve.py) -----------------------

def _normalize_image(img: np.ndarray) -> np.ndarray:
    p1, p99 = np.percentile(img, (1.0, 99.0))
    if p99 <= p1:
        return img.astype(np.float32) * 0
    img = img.clip(p1, p99)
    return ((img - p1) / (p99 - p1)).astype(np.float32)


def _astroalign_worker(src_gray: np.ndarray, fits_norm: np.ndarray, conn: Any) -> None:
    import astroalign as aa

    try:
        transform, (src_pts, _) = aa.find_transform(
            src_gray, fits_norm, max_control_points=100, detection_sigma=5, min_area=5
        )
        conn.send(("ok", transform.params, float(transform.scale), float(transform.rotation), len(src_pts)))
    except Exception as exc:
        conn.send(("error", str(exc)))
    finally:
        conn.close()


def _find_transform(plate_img: Image.Image, fits_norm: np.ndarray) -> dict[str, Any] | None:
    """Match stars between the photo and the FITS, giving up after ASTROALIGN_TIMEOUT_S.

    Runs in a child process so a hopeless attempt (usually the wrong mirror
    orientation) can be killed. Lambda has no shared memory, so this uses a Pipe,
    not a Queue; fork (Linux) lets the child reuse already-loaded modules.
    """
    # Photos have dark stars on a light background; invert to match the FITS.
    src_gray = 1.0 - _normalize_image(np.array(plate_img, dtype=np.float32))

    method = "fork" if "fork" in mp.get_all_start_methods() else "spawn"  # spawn: Windows testing
    ctx = mp.get_context(method)
    recv_end, send_end = ctx.Pipe(duplex=False)
    proc = ctx.Process(target=_astroalign_worker, args=(src_gray, fits_norm, send_end))
    proc.start()
    send_end.close()

    got_result = recv_end.poll(ASTROALIGN_TIMEOUT_S)
    if not got_result:
        LOG.info("Astroalign timed out after %.1fs", ASTROALIGN_TIMEOUT_S)
        proc.terminate()
        proc.join()
        return None
    try:
        status, *payload = recv_end.recv()
    except EOFError:
        status, payload = "error", ["worker exited without a result"]
    proc.join()

    if status == "error":
        LOG.info("Astroalign failed: %s", payload[0])
        return None
    params, scale, rotation, n = payload
    LOG.info("Astroalign: matches=%d scale=%.5f rot_deg=%.3f", n, scale, np.degrees(rotation))
    return {"params": np.asarray(params), "scale": scale, "rotation": rotation, "matches": n}


def _cover_resize(img: Image.Image, target_w: int, target_h: int) -> tuple[Image.Image, float]:
    """Scale img so it covers (target_w, target_h); return (resized, scale)."""
    w, h = img.size
    scale = max(target_w / max(1, w), target_h / max(1, h))
    resized = img.resize(
        (max(1, int(round(w * scale))), max(1, int(round(h * scale)))),
        resample=Image.Resampling.BICUBIC,
    )
    return resized, float(scale)


def _compose_wcs_from_affine(reference_wcs: WCS, affine_a: np.ndarray, affine_t: np.ndarray) -> dict[str, Any]:
    """Turn a photo->FITS pixel transform plus the FITS WCS into a WCS header for the photo."""
    if reference_wcs.wcs.has_cd():
        ref_cd = np.array(reference_wcs.wcs.cd, dtype=np.float64)
    else:
        ref_cd = np.array(reference_wcs.wcs.get_pc(), dtype=np.float64) @ np.diag(reference_wcs.wcs.cdelt)

    ones = np.ones(2)
    cd_src = ref_cd @ affine_a
    b = affine_t + ones - affine_a @ ones
    crpix_src = np.linalg.solve(affine_a, np.array(reference_wcs.wcs.crpix, dtype=np.float64) - b)

    header = reference_wcs.to_header(relax=True)
    for key in list(header.keys()):
        if key.startswith(("PC", "CD")) or key in ("CDELT1", "CDELT2"):
            del header[key]
    header["CRPIX1"], header["CRPIX2"] = float(crpix_src[0]), float(crpix_src[1])
    header["CD1_1"], header["CD1_2"] = float(cd_src[0, 0]), float(cd_src[0, 1])
    header["CD2_1"], header["CD2_2"] = float(cd_src[1, 0]), float(cd_src[1, 1])
    return dict(header)


def solve_plate(plate_id: str) -> dict[str, Any]:
    """Align the plate JPG to its FITS mosaic and return the JPG's WCS."""
    photo = read_plate_photo(plate_id)
    reference_wcs, fits_data = read_plate_fits(plate_id)
    target_h, target_w = fits_data.shape
    fits_norm = _normalize_image(fits_data)
    src_w, src_h = photo.size

    # Resize/crop the photo to roughly the FITS size so star patterns line up.
    resized, resize_scale = _cover_resize(photo, target_w, target_h)
    x0 = max(0, min((resized.width - target_w) // 2, resized.width - target_w))
    y0 = max(0, min((resized.height - target_h) // 2, resized.height - target_h))
    match_img = resized.crop((x0, y0, x0 + target_w, y0 + target_h))

    # Mirrored first: most plate photos are mirror images of the scans.
    candidates = [
        ("source_mirror", match_img.transpose(Image.Transpose.FLIP_LEFT_RIGHT)),
        ("source", match_img),
    ]
    for label, candidate in candidates:
        t = _find_transform(candidate, fits_norm)
        if t is None:
            continue

        # Pixel map from original photo -> the resized/cropped (and maybe mirrored) image.
        if label == "source_mirror":
            base_a = np.array([[-1.0, 0.0], [0.0, 1.0]]) * resize_scale
            base_t = np.array([target_w - 1.0 + x0, -y0])
        else:
            base_a = np.eye(2) * resize_scale
            base_t = np.array([-x0, -y0], dtype=np.float64)

        a, tt = t["params"][:2, :2], t["params"][:2, 2]
        wcs_header = _compose_wcs_from_affine(reference_wcs, a @ base_a, tt + a @ base_t)
        return {
            "plate_id": plate_id,
            "wcs_header": wcs_header,
            "image_width_px": src_w,
            "image_height_px": src_h,
            "alignment_meta": {
                "backend": "astroalign",
                "parity": label,
                "matches": t["matches"],
                "transform_scale": t["scale"],
                "transform_rotation_deg": float(np.degrees(t["rotation"])),
            },
        }

    raise RuntimeError("astroalign found no star match in either orientation")


# --- WWT link and WTML (WTML structure from src/dasch_sky_mosaic/wtml_gen.py) ---

def _own_url(event: dict[str, Any]) -> str | None:
    """This function's public URL, from the request (absent in console test events)."""
    domain = (event.get("requestContext") or {}).get("domainName")
    return f"https://{domain}/" if domain else None


def build_wwt_link(event: dict[str, Any], plate_id: str) -> str | None:
    """WWT webclient link that loads this plate's WTML from this function."""
    own_url = _own_url(event)
    if own_url is None:
        return None
    wtml_url = f"{own_url}?{urlencode({'plate': plate_id, 'format': 'wtml'})}"
    return f"{WEBCLIENT_URL}?{urlencode({'wtml': wtml_url})}"


def build_wtml(plate_id: str, result: dict[str, Any]) -> str:
    """WTML (a Folder with one Place) showing the plate JPG at its solved position."""
    header = result["wcs_header"]
    width, height = int(result["image_width_px"]), int(result["image_height_px"])

    # wwt_data_formats does the WCS -> WWT placement maths (scale, rotation, parity).
    placement = WwtImageSet()
    placement.tile_levels = 0
    placement.set_position_from_wcs(header, width=width, height=height)

    deg_per_px = float(placement.base_degrees_per_tile)
    center_ra = float(placement.center_x)
    center_dec = float(placement.center_y)

    image_url = s3.generate_presigned_url(
        "get_object",
        Params={
            "Bucket": PLATES_BUCKET,
            "Key": _photo_key(plate_id),
            # The plate JPGs are stored as binary/octet-stream; tell S3 to serve
            # this one as an image to display, not a file to download.
            "ResponseContentType": "image/jpeg",
            "ResponseContentDisposition": "inline",
        },
        ExpiresIn=IMAGE_URL_EXPIRY_S,
    )

    folder = ET.Element(
        "Folder", {"Name": f"DASCH {plate_id}", "Group": "Explorer", "Type": "Sky", "Searchable": "True"}
    )
    place = ET.SubElement(
        folder,
        "Place",
        {
            "Name": plate_id,
            "DataSetType": "Sky",
            "RA": f"{center_ra / 15.0:.8f}",  # hours
            "Dec": f"{center_dec:.8f}",
            "ZoomLevel": f"{max(width, height) * deg_per_px * 6.0:.6f}",
            "Opacity": "100",
        },
    )
    foreground = ET.SubElement(place, "ForegroundImageSet")
    imageset = ET.SubElement(
        foreground,
        "ImageSet",
        {
            "Generic": "False",
            "DataSetType": "Sky",
            "BandPass": "Visible",
            "Name": f"DASCH {plate_id}",
            "Url": image_url,  # ElementTree escapes the & characters in the XML
            "TileLevels": "0",
            "BaseTileLevel": "0",
            "BaseDegreesPerTile": f"{deg_per_px:.12f}",
            "FileType": ".jpg",
            "Projection": "SkyImage",
            "BottomsUp": "True" if bool(placement.bottoms_up) else "False",
            "CenterX": f"{center_ra:.8f}",
            "CenterY": f"{center_dec:.8f}",
            "OffsetX": f"{float(placement.offset_x):.8f}",
            "OffsetY": f"{float(placement.offset_y):.8f}",
            "Rotation": f"{float(placement.rotation_deg):.8f}",
            "Sparse": "False",
            "WidthFactor": f"{float(placement.width_factor):.8f}",
            "QuadTreeMap": "",
        },
    )
    ET.SubElement(imageset, "ThumbnailUrl")
    ET.SubElement(imageset, "Credits").text = CREDITS
    ET.SubElement(imageset, "CreditsUrl").text = CREDITS_URL
    ET.SubElement(imageset, "Description").text = f"DASCH plate {plate_id}"

    return "<?xml version='1.0' encoding='UTF-8'?>\n" + ET.tostring(folder, encoding="unicode")


