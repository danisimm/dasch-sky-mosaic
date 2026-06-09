import json
from pathlib import Path
from urllib.parse import quote
from astropy.io import fits

out = Path('data/output/wwt_browser_links.json')
links = json.loads(out.read_text(encoding='utf-8'))

for plate in ('ai36023','ai42221'):
    j = Path(f'data/output/{plate}_wtml.json')
    if not j.exists():
        print('missing', j)
        continue
    data = json.loads(j.read_text(encoding='utf-8'))
    p = data['paired'][0]
    image = Path(p['aligned_photo_path'])
    image_uri = image.as_uri()
    mosaic = Path(p['wcs_mosaic_path'])
    # ensure fits (funpacked) path
    if mosaic.suffix == '.fz':
        fits_path = mosaic.with_suffix('')
        if fits_path.suffix == '.fit':
            fits_path = fits_path.with_suffix('.fits')
    else:
        fits_path = mosaic
    with fits.open(fits_path) as hdul:
        header = hdul[0].header
        ra = header.get('CRVAL1') or header.get('CRVAL') or 0.0
        dec = header.get('CRVAL2') or 0.0
    placement = p['skyimage_placement']
    scale = max(1e-6, placement['base_degrees_per_tile'])
    x = placement['offset_x']
    y = placement['offset_y']
    rotation = placement['rotation_deg'] % 360.0
    show = (
        'http://www.worldwidetelescope.org/wwtweb/ShowImage.aspx?reverseparity=False'
        f'&scale={scale:.5f}&name={plate}&imageurl={quote(image_uri, safe="")}'
        '&credits=DASCH+%2F+Harvard+College+Observatory&creditsUrl=https%3A%2F%2Fdasch.cfa.harvard.edu%2F'
        f'&ra={ra:.9f}&dec={dec:.9f}&y={y:.8f}&x={x:.8f}&rotation={rotation:.8f}&wtml=true'
    )
    launch = 'http://www.worldwidetelescope.org/webclient/default.aspx?wtml=' + quote(show, safe='')
    links.append({'plate_id': plate, 'showimage_url': show, 'launch_url': launch})

out.write_text(json.dumps(links, indent=2), encoding='utf-8')
print('updated', out)
