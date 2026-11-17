"""Read the cfg1000 handoff and cfg1020 diagonal queue without modifying it."""
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
import json
import shutil

BASE = Path('/Users/admin/Downloads/MetalClover-results')
WORK = Path('/private/tmp/chroma-cfg1020-production')


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        return {}


jobs = []
for configuration, times in [(1000, [9, 10]), (1020, range(1, 11))]:
    for t in times:
        point = BASE / f'cfg{configuration}-t{t}-unsmeared'
        queue = read(str(point) + '.queue.json')
        pipeline = read(str(point) + '.pipeline.json')
        run = read(str(point) + '.run.json')
        jobs.append({'configuration': configuration, 'source_sink_time': t,
                     'status': queue.get('status', 'not_started'),
                     'stage': pipeline.get('stage', 'waiting'),
                     'completed_RHS': run.get('completed_rhs', 0), 'total_RHS': 20736,
                     'maximum_full_relative': run.get('max_full_relative'),
                     'waiting_for': queue.get('waiting_for')})
print(json.dumps({'local_time': datetime.now(ZoneInfo('America/New_York')).isoformat(),
                  'free_disk_GiB': shutil.disk_usage(BASE).free / 2**30,
                  'storage': read(WORK / 'storage-supervisor.json'),
                  'jobs': jobs}, indent=2))
