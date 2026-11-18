import numpy as np
import sep
sep.set_extract_pixstack(1000000)
from astropy.io import fits
from scipy.spatial import cKDTree

FITS_PATH = 'data/cache/scamp_solve/b02312/plate_rot0.fits'
LDAC_PATH = 'data/cache/scamp_solve/b02312/plate_rot0.ldac'

with fits.open(FITS_PATH) as h:
    img = h[0].data.astype('float32')

bkg = sep.Background(img)
thresh = 5 * bkg.globalrms
aa_sources = sep.extract(img - bkg.back(), thresh, minarea=5)
aa_xy = np.column_stack([aa_sources['x'], aa_sources['y']])
aa_flux = aa_sources['flux']
print(f'astroalign-style (sep) detections: {len(aa_xy)}')

with fits.open(LDAC_PATH) as h:
    obj = h['LDAC_OBJECTS'].data
sex_xy = np.column_stack([obj['XWIN_IMAGE'] - 1, obj['YWIN_IMAGE'] - 1])  # FITS 1-indexed -> 0-indexed
sex_flux = obj['FLUX_AUTO']
print(f'SExtractor detections: {len(sex_xy)}')

tol_px = 3.0
tree = cKDTree(aa_xy)
dist, idx = tree.query(sex_xy, k=1)
matched = dist < tol_px
print(f'SExtractor sources with an astroalign-style match within {tol_px}px: {matched.sum()} / {len(sex_xy)} ({100*matched.mean():.1f}%)')

tree2 = cKDTree(sex_xy)
dist2, idx2 = tree2.query(aa_xy, k=1)
matched2 = dist2 < tol_px
print(f'astroalign-style sources with a SExtractor match within {tol_px}px: {matched2.sum()} / {len(aa_xy)} ({100*matched2.mean():.1f}%)')

print()
print('Top 10 brightest astroalign-style sources (x,y,flux):')
order = np.argsort(aa_flux)[::-1][:10]
for i in order:
    print(f'  ({aa_xy[i,0]:.1f}, {aa_xy[i,1]:.1f})  flux={aa_flux[i]:.1f}')

print()
print('Top 10 brightest SExtractor sources (x,y,flux):')
order2 = np.argsort(sex_flux)[::-1][:10]
for i in order2:
    print(f'  ({sex_xy[i,0]:.1f}, {sex_xy[i,1]:.1f})  flux={sex_flux[i]:.1f}')
