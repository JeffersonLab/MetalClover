"""Read the local queue; never start or modify a numerical job."""
from pathlib import Path
import datetime
import json
import shutil

BASE = Path('/Users/admin/Downloads/MetalClover-results')


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        return {}


jobs = []
for t in range(4, 21):
    prefix = BASE / f'cfg1000-t{t}-unsmeared'
    queue = read(str(prefix) + '.queue.json')
    pipeline = read(str(prefix) + '.pipeline.json')
    run = read(str(prefix) + '.run.json')
    jobs.append({'source_sink_time': t,
                 'queue_status': queue.get('status', 'not_started'),
                 'stage': pipeline.get('stage', queue.get('stage', 'waiting')),
                 'queue_PID': queue.get('queue_PID'),
                 'pipeline_PID': queue.get('pipeline_PID'),
                 'completed_rhs': run.get('completed_rhs', 0),
                 'total_rhs': 20736,
                 'waiting_for': queue.get('waiting_for')})
print(json.dumps({'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  'free_disk_GiB': shutil.disk_usage(BASE).free / 2**30,
                  'jobs': jobs}, indent=2))
