from pathlib import Path
import sys
from dasch_sky_mosaic.fetch import download_mosaic_paths

if len(sys.argv) < 2:
    print('Usage: download_mosaic.py <plate_id>...')
    raise SystemExit(2)
plates = [p.lower() for p in sys.argv[1:]]
paths = download_mosaic_paths(Path('data/cache/dasch_session'), plates, binning=16, api_base='auto', api_key=None)
for k,v in paths.items():
    print(k, v)
