import sys
sys.path.insert(0, 'src')
import json
import re
import numpy as np
import pandas as pd
from pathlib import Path
from astropy.coordinates import SkyCoord
import astropy.units as u

from dasch_sky_mosaic.wcs_solve import solve_plate_wcs

exp = pd.read_csv('dev/exposures.csv')

photo_dir = Path('data/cache/photos')
fits_dir = Path('data/cache/mosaic_package')
wcs_out_dir = Path('data/output/wcs')
wcs_out_dir.mkdir(parents=True, exist_ok=True)

plate_ids = sorted(p.name.split('_pphoto_all')[0] for p in photo_dir.glob('*.jpg'))
fits_by_plate = {}
for f in fits_dir.glob('*.fits'):
    m = re.match(r'(.+?)_mosaic_(\d+)_16', f.name)
    if m:
        fits_by_plate.setdefault(m.group(1), []).append((int(m.group(2)), f))

def pa_of(CD):
    v_y = CD @ np.array([0.0, 1.0])
    return np.degrees(np.arctan2(v_y[0], v_y[1]))

results = []
for plate_id in plate_ids:
    fits_candidates = fits_by_plate.get(plate_id)
    if not fits_candidates:
        print(f'{plate_id}: no fits cached, skip')
        continue
    solnum, fits_path = sorted(fits_candidates)[0]
    photo_path = photo_dir / f'{plate_id}_pphoto_all.jpg'
    try:
        res = solve_plate_wcs(plate_id, photo_path, fits_path)
    except Exception as exc:
        print(f'{plate_id}: astroalign FAILED - {exc}')
        results.append(dict(plate_id=plate_id, solved=False, error=str(exc)))
        continue

    (wcs_out_dir / f'{plate_id}_wcs.json').write_text(json.dumps(res, indent=2))

    rows = exp[(exp['plateid'] == plate_id) & (exp['deleted'] != 1.0)]
    pref = rows[(rows['solutionNumber'] == solnum) & (rows['exposureNumber'] == 0)]
    if pref.empty:
        pref = rows[rows['exposureNumber'] == 0]
    if pref.empty:
        pref = rows
    if pref.empty:
        print(f'{plate_id}: no exposures.csv row, skip')
        results.append(dict(plate_id=plate_id, solved=True, error='no_exposures_row'))
        continue
    row = pref.iloc[0]

    gt = res['wcs_header']
    CD_actual = np.array([[gt['CD1_1'], gt['CD1_2']], [gt['CD2_1'], gt['CD2_2']]])
    CD_raw = np.array([[row['delta_ra_x'], row['delta_ra_y']], [row['delta_dec_x'], row['delta_dec_y']]])

    jpg_w, jpg_h = res['image_width_px'], res['image_height_px']
    naxis1, naxis2 = row['naxis1'], row['naxis2']

    guess_scale = np.sqrt(abs(np.linalg.det(CD_raw)) * naxis1 * naxis2 / (jpg_w * jpg_h))
    actual_scale = np.sqrt(abs(np.linalg.det(CD_actual)))

    pa_actual = pa_of(CD_actual)
    pa_raw = pa_of(CD_raw)
    parity_actual = np.sign(np.linalg.det(CD_actual))
    parity_raw = np.sign(np.linalg.det(CD_raw))

    c1 = SkyCoord(gt['CRVAL1'] * u.deg, gt['CRVAL2'] * u.deg)
    c2 = SkyCoord(row['ctr_ra'] * u.deg, row['ctr_dec'] * u.deg)
    sep_arcmin = c1.separation(c2).arcmin

    results.append(dict(
        plate_id=plate_id,
        solved=True,
        error=None,
        solnum_used=solnum,
        exp_solnum=int(row['solutionNumber']),
        sep_arcmin=sep_arcmin,
        guess_scale_degpx=guess_scale,
        actual_scale_degpx=actual_scale,
        scale_ratio=actual_scale / guess_scale,
        pa_actual_deg=pa_actual,
        pa_raw_deg=pa_raw,
        pa_diff_deg=(pa_actual - pa_raw) % 360,
        pa_diff_mod90=(pa_actual - pa_raw) % 90,
        parity_match=bool(parity_actual == parity_raw),
        n_matches=res['alignment_meta']['matches'],
        naxis1=naxis1, naxis2=naxis2, jpg_w=jpg_w, jpg_h=jpg_h,
    ))

df = pd.DataFrame(results)
out_csv = Path('data/output/analysis/scamp_seed_analysis.csv')
df.to_csv(out_csv, index=False)
pd.set_option('display.width', 220)
pd.set_option('display.max_columns', 30)
print(df.to_string(index=False))

solved = df[df['solved'] == True].dropna(subset=['sep_arcmin'])
print()
print(f'n solved with exposures.csv match: {len(solved)} / {len(df)}')
print('sep_arcmin: mean=%.2f median=%.2f max=%.2f' % (solved.sep_arcmin.mean(), solved.sep_arcmin.median(), solved.sep_arcmin.max()))
print('scale_ratio: mean=%.3f std=%.3f min=%.3f max=%.3f' % (solved.scale_ratio.mean(), solved.scale_ratio.std(), solved.scale_ratio.min(), solved.scale_ratio.max()))
print()
print('pa_diff_mod90 (residual after removing nearest-90-deg step):')
print(solved[['plate_id', 'pa_diff_deg', 'pa_diff_mod90', 'parity_match']].to_string(index=False))
print()
print('Saved analysis to', out_csv)
