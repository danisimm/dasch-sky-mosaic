from __future__ import annotations

import io
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any
from urllib.parse import urlencode
from urllib.request import Request, urlopen

import requests
from astropy.io import fits
from astropy.wcs import WCS


import logging


LOG = logging.getLogger(__name__)


@dataclass(frozen=True)
class AstrometrySolveConfig:
    api_key: str
    api_base: str
    working_dir: Path
    timeout_seconds: float
    upload_timeout_seconds: float
    poll_interval_seconds: float
    overwrite: bool
    ra_hint_deg: float | None = None
    dec_hint_deg: float | None = None
    radius_hint_deg: float | None = None
    scale_low_arcsec_per_pix: float | None = None
    scale_high_arcsec_per_pix: float | None = None


def _read_wcs_header(wcs_path: Path) -> dict[str, Any]:
    with fits.open(wcs_path, memmap=True) as hdul:
        for hdu in hdul:
            header = hdu.header
            # A valid astrometric WCS should provide celestial axis types.
            ctype1 = str(header.get("CTYPE1", ""))
            ctype2 = str(header.get("CTYPE2", ""))
            if "RA" in ctype1 and "DEC" in ctype2:
                wcs = WCS(header)
                return dict(wcs.to_header(relax=True))
    raise RuntimeError(f"astrometry.net did not produce a valid celestial WCS in {wcs_path}")


def _post_request_json(
    http: requests.Session,
    api_base: str,
    endpoint: str,
    payload: dict[str, Any],
    timeout_seconds: float,
    files: dict[str, Any] | None = None,
) -> dict[str, Any]:
    url = f"{api_base.rstrip('/')}/{endpoint.lstrip('/')}"
    response = http.post(
        url,
        data={"request-json": json.dumps(payload)},
        files=files,
        timeout=timeout_seconds,
    )
    response.raise_for_status()
    result = response.json()
    if not isinstance(result, dict):
        raise RuntimeError(f"invalid API response from {url}: expected object")
    return result


def _post_json(
    http: requests.Session,
    url: str,
    timeout_seconds: float,
    payload: dict[str, Any] | None = None,
) -> dict[str, Any]:
    request = Request(
        url,
        data=urlencode({"request-json": json.dumps(payload or {})}).encode("utf-8"),
        headers={"Content-Type": "application/x-www-form-urlencoded", "Accept": "application/json"},
    )
    with urlopen(request, timeout=timeout_seconds) as response:
        raw = response.read().decode("utf-8")
    result = json.loads(raw)
    if not isinstance(result, dict):
        raise RuntimeError(f"invalid API response from {url}: expected object")
    return result


def _api_root_from_base(api_base: str) -> str:
    base = api_base.rstrip("/")
    if base.endswith("/api"):
        return base[: -len("/api")]
    return base


def solve_image_wcs(
    image_path: Path,
    plate_id: str,
    config: AstrometrySolveConfig,
) -> tuple[dict[str, Any], dict[str, Any]]:
    """Run astrometry.net Nova API solve and return a FITS-compatible WCS header.

    Returns (wcs_header, solver_meta).
    """
    config.working_dir.mkdir(parents=True, exist_ok=True)

    stem = f"{plate_id.lower()}_solve"
    wcs_path = config.working_dir / f"{stem}.wcs.fits"
    submit_meta_path = config.working_dir / f"{stem}.submission.json"

    if wcs_path.exists() and not config.overwrite:
        return _read_wcs_header(wcs_path), {
            "status": "cached",
            "wcs_path": str(wcs_path),
            "submission_meta": str(submit_meta_path),
            "job_id": None,
            "sub_id": None,
        }

    api_base = config.api_base.rstrip("/")
    api_root = _api_root_from_base(api_base)

    with requests.Session() as http:
        LOG.info("Astrometry: logging in to %s", api_root)
        login = _post_request_json(
            http=http,
            api_base=api_base,
            endpoint="login",
            payload={"apikey": config.api_key},
            timeout_seconds=config.upload_timeout_seconds,
        )
        if login.get("status") != "success" or not login.get("session"):
            raise RuntimeError(f"astrometry.net login failed: {login}")
        session_key = str(login["session"])
        LOG.info("Astrometry: login successful")

        upload_payload: dict[str, Any] = {
            "session": session_key,
            "allow_commercial_use": "d",
            "allow_modifications": "d",
            "publicly_visible": "n",
        }
        if config.ra_hint_deg is not None and config.dec_hint_deg is not None:
            upload_payload["center_ra"] = float(config.ra_hint_deg)
            upload_payload["center_dec"] = float(config.dec_hint_deg)
            if config.radius_hint_deg is not None:
                upload_payload["radius"] = float(config.radius_hint_deg)
        if config.scale_low_arcsec_per_pix is not None and config.scale_high_arcsec_per_pix is not None:
            upload_payload["scale_units"] = "arcsecperpix"
            upload_payload["scale_type"] = "ul"
            upload_payload["scale_lower"] = float(config.scale_low_arcsec_per_pix)
            upload_payload["scale_upper"] = float(config.scale_high_arcsec_per_pix)

        with image_path.open("rb") as f_in:
            LOG.info("Astrometry: uploading %s", image_path.name)
            upload = _post_request_json(
                http=http,
                api_base=api_base,
                endpoint="upload",
                payload=upload_payload,
                timeout_seconds=config.upload_timeout_seconds,
                files={
                    "file": (
                        image_path.name,
                        f_in,
                        "application/octet-stream",
                    )
                },
            )

        if upload.get("status") != "success" or upload.get("subid") is None:
            raise RuntimeError(f"astrometry.net upload failed: {upload}")
        sub_id = int(upload["subid"])
        LOG.info("Astrometry: submission id %d", sub_id)

        deadline = time.monotonic() + config.timeout_seconds
        job_id: int | None = None
        while time.monotonic() < deadline:
            LOG.info("Astrometry: polling submission %d", sub_id)
            sub_status = _post_json(
                http,
                f"{api_base}/submissions/{sub_id}",
                timeout_seconds=config.upload_timeout_seconds,
            )
            jobs = sub_status.get("jobs")
            if isinstance(jobs, list):
                valid_jobs = [int(j) for j in jobs if isinstance(j, int) or (isinstance(j, str) and str(j).isdigit())]
                if valid_jobs:
                    job_id = valid_jobs[-1]
                    LOG.info("Astrometry: selected job id %d", job_id)
                    break
            time.sleep(config.poll_interval_seconds)

        if job_id is None:
            raise RuntimeError(f"astrometry.net timed out waiting for job creation (subid={sub_id})")

        final_job_status: dict[str, Any] | None = None
        while time.monotonic() < deadline:
            LOG.info("Astrometry: polling job %d", job_id)
            job_status = _post_json(
                http,
                f"{api_base}/jobs/{job_id}",
                timeout_seconds=config.upload_timeout_seconds,
            )
            status = str(job_status.get("status", "")).lower()
            if status == "success":
                final_job_status = job_status
                break
            if status == "failure":
                raise RuntimeError(f"astrometry.net solve failed for job {job_id}")
            time.sleep(config.poll_interval_seconds)

        if final_job_status is None:
            raise RuntimeError(f"astrometry.net timed out waiting for solve completion (job={job_id})")

        wcs_url = f"{api_root}/wcs_file/{job_id}"
        LOG.info("Astrometry: downloading WCS from %s", wcs_url)
        wcs_response = http.get(wcs_url, timeout=config.upload_timeout_seconds)
        wcs_response.raise_for_status()
        wcs_bytes = wcs_response.content
        if not wcs_bytes:
            raise RuntimeError(f"astrometry.net returned empty WCS file for job {job_id}")
        wcs_path.write_bytes(wcs_bytes)

        calibration: dict[str, Any] | None = None
        try:
            calibration = _post_json(
                http,
                f"{api_base}/jobs/{job_id}/calibration",
                timeout_seconds=config.upload_timeout_seconds,
            )
        except Exception:
            calibration = None

    wcs_header = _read_wcs_header(wcs_path)
    submit_meta = {
        "status": "solved",
        "sub_id": sub_id,
        "job_id": job_id,
        "wcs_url": f"{api_root}/wcs_file/{job_id}",
        "wcs_path": str(wcs_path),
        "calibration": calibration,
    }
    submit_meta_path.write_text(json.dumps(submit_meta, indent=2), encoding="utf-8")
    return wcs_header, submit_meta
