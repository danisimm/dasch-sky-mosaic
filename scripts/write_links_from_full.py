import json
from pathlib import Path
from urllib.parse import quote

out = Path('data/output/wwt_browser_links.json')
links = json.loads(out.read_text(encoding='utf-8'))

def make_entry(plate):
    j = Path(f'data/output/{plate}_wtml_full.json')
    data = json.loads(j.read_text(encoding='utf-8'))
    p = data['paired'][0]
    image_uri = Path(p['aligned_photo_path']).as_uri()
    pl = p['skyimage_placement']
    ra = p['center_ra_deg']
    dec = p['center_dec_deg']
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
    return {'plate_id': plate, 'showimage_url': show, 'launch_url': launch}

for plate in ('ai36023','ai42221'):
    links = [l for l in links if l.get('plate_id') != plate]
    links.append(make_entry(plate))

out.write_text(json.dumps(links, indent=2), encoding='utf-8')
print('written', out)
