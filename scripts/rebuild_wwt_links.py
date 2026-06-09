import json
from pathlib import Path
from urllib.parse import quote

plates = ['a11740','br00130','ax01822','ai36023','ai42221']
out = Path('data/output/wwt_browser_links.json')
links = []

for plate in plates:
    j = Path(f'data/output/{plate}_wtml_full.json')
    if not j.exists():
        j = Path(f'data/output/{plate}_wtml.json')
    if not j.exists():
        print('missing', plate)
        continue
    data = json.loads(j.read_text(encoding='utf-8'))
    p = data['paired'][0]
    image_uri = Path(p['aligned_photo_path']).as_uri()
    pl = p['skyimage_placement']
    ra = p.get('center_ra_deg') or p.get('skyimage_placement',{}).get('center_x')
    dec = p.get('center_dec_deg') or p.get('skyimage_placement',{}).get('center_y')
    x = pl['offset_x']
    y = pl['offset_y']
    rotation = pl['rotation_deg'] % 360.0
    scale = pl['base_degrees_per_tile']
    show = (
        'http://www.worldwidetelescope.org/wwtweb/ShowImage.aspx?reverseparity=False'
        f'&scale={scale:.5f}&name={plate}&imageurl={quote(image_uri, safe="")}'
        '&credits=DASCH+%2F+Harvard+College+Observatory&creditsUrl=https%3A%2F%2Fdasch.cfa.harvard.edu%2F'
        f'&ra={ra:.9f}&dec={dec:.9f}&y={y:.8f}&x={x:.8f}&rotation={rotation:.8f}&wtml=true'
    )
    launch = 'http://www.worldwidetelescope.org/webclient/default.aspx?wtml=' + quote(show, safe='')
    links.append({'plate_id': plate, 'showimage_url': show, 'launch_url': launch})

out.write_text(json.dumps(links, indent=2), encoding='utf-8')
print('rebuilt', out)
