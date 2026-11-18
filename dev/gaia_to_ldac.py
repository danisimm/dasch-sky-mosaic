"""One-off: convert the local Gaia catalog into SCAMP's reference-catalog format.

SCAMP's -ASTREFCAT_NAME (ASTREF_CATALOG FILE mode) expects a FITS binary
table with specific recognized column names (X_WORLD/Y_WORLD for position,
ERRA_WORLD/ERRB_WORLD for position error, MAG for magnitude) - best
understanding from documentation, not yet verified against a real SCAMP run,
so this may need adjusting once tested.

Proper motions are dropped here rather than propagated: DASCH plates span
~1880s-1990s with no single representative epoch, and FILE-mode reference
catalogs don't have a documented way to carry per-source epoch info for SCAMP
to propagate itself. High-proper-motion stars may simply fail to cross-match
at the current epoch and drop out of the fit - acceptable given the G<14 cut
already has density margin (see the magnitude cutoff discussion).

Run once: python dev/gaia_to_ldac.py
"""
import numpy as np
from astropy.io import fits
from astropy.table import Table

SRC = "data/reference_catalogs/gaia_dr3_g14.fits"
DEST = "data/reference_catalogs/gaia_dr3_g14_ldac.fits"

# Gaia DR3 astrometric precision at these magnitudes is far better than
# anything else in this pipeline (our own initial-guess and plate-detector
# errors dominate) - a small constant stand-in is fine, SCAMP just needs a
# nonzero error to weight by.
ASSUMED_POSITION_ERROR_DEG = 1e-6


def main() -> None:
    table = Table.read(SRC)

    columns = [
        fits.Column(name="X_WORLD", format="D", unit="deg", array=table["ra"]),
        fits.Column(name="Y_WORLD", format="D", unit="deg", array=table["dec"]),
        fits.Column(name="ERRA_WORLD", format="D", unit="deg", array=np.full(len(table), ASSUMED_POSITION_ERROR_DEG)),
        fits.Column(name="ERRB_WORLD", format="D", unit="deg", array=np.full(len(table), ASSUMED_POSITION_ERROR_DEG)),
        fits.Column(name="MAG", format="E", array=table["phot_g_mean_mag"]),
    ]
    objects_hdu = fits.BinTableHDU.from_columns(columns, name="LDAC_OBJECTS")
    fits.HDUList([fits.PrimaryHDU(), objects_hdu]).writeto(DEST, overwrite=True)
    print(f"wrote {len(table)} sources -> {DEST}")


if __name__ == "__main__":
    main()
