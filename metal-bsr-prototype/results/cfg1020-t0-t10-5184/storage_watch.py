"""Maintain a full-sequence disk gate without deleting or moving any data."""
from pathlib import Path
import argparse
import datetime
import json
import os
import shutil
import signal
import time

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
REPORT = ROOT / 'metal-bsr-prototype/results/cfg1020-t0-t10-5184'
WORK = Path('/private/tmp/chroma-cfg1020-production')
BASE = Path('/Users/admin/Downloads/MetalClover-results')
GIB = 2**30
POINT = 6879707264
STANDARD = 16 * 5184**2 * 16 + 64 * 2**20
EVIDENCE = GIB // 4
PROJECTION_REQUIRED_GIB = 1.1 * (16 * 5184**2 * 16) / GIB + 4


def read(path):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        return {}


def write(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp-{os.getpid()}')
    with temporary.open('w') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def allocation(path):
    try:
        info = path.stat()
        assert not path.is_symlink()
        return [info.st_ino, info.st_size, info.st_blocks * 512]
    except FileNotFoundError:
        return [None, 0, 0]


def snapshot():
    states = []
    for cfg, t in [(1000, 9), (1000, 10)] + [(1020, t) for t in range(1, 11)]:
        point = BASE / f'cfg{cfg}-t{t}-unsmeared'
        standard = BASE / f'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb{cfg}'
        queue = read(str(point) + '.queue.json')
        states.append({'configuration': cfg, 'time': t, 'queue_status': queue.get('status', 'not_started'),
                       'pipeline_exists': Path(str(point) + '.pipeline.json').exists(),
                       'point': allocation(Path(str(point) + '.npy')),
                       'standard': allocation(standard)})
    return states


def forecast(free_bytes, states):
    free = free_bytes - GIB
    events = [{'event': 'reserve 1 GiB for unrelated modest disk use', 'free_GiB': free / GIB}]
    starts, projections = [], []
    for state in states:
        if state['queue_status'] == 'complete':
            continue
        t = state['time']
        label = f"cfg{state.get('configuration', 1020)} time{t}"
        if not state['pipeline_exists']:
            starts.append(free / GIB)
        events.append({'event': f'{label} remaining work begins', 'free_GiB': free / GIB})
        # Existing allocation is already charged to the measured free space.
        point_allocation = max(POINT, state['point'][2])
        free -= max(0, POINT - state['point'][2])
        free -= EVIDENCE
        if state['standard'][0] is None:
            projections.append(free / GIB)
        free -= max(0, STANDARD - state['standard'][2])
        events.append({'event': f'{label} peak output before validated point retirement',
                       'free_GiB': free / GIB})
        free += point_allocation
        events.append({'event': f'{label} accepted point retired; standard retained',
                       'free_GiB': free / GIB})
    minimum = min(event['free_GiB'] for event in events)
    minimum_start = min(starts, default=18.0)
    minimum_projection = min(projections, default=PROJECTION_REQUIRED_GIB)
    deficit = max(0, 4 - minimum, 18 - minimum_start,
                  PROJECTION_REQUIRED_GIB - minimum_projection)
    fits = deficit == 0
    return {'status': 'pass' if fits else 'insufficient_space',
            'forecast_all_runs_fit': fits,
            'all_standard_perambulators_retained': True,
            'actual_free_disk_GiB': free_bytes / GIB,
            'additional_free_space_required_GiB': deficit,
            'required_actual_free_disk_GiB': free_bytes / GIB + deficit,
            'minimum_forecast_free_GiB': minimum,
            'minimum_new_job_start_free_GiB': minimum_start,
            'minimum_before_projection_GiB': minimum_projection,
            'native_projection_required_free_GiB': PROJECTION_REQUIRED_GIB,
            'final_forecast_free_GiB': free / GIB,
            'standard_budget_bytes_each': STANDARD,
            'log_evidence_reserve_bytes_per_remaining_time': EVIDENCE,
            'global_reserve_bytes': GIB,
            'no_files_deleted_or_moved_by_storage_supervisor': True,
            'events': events, 'snapshot': states}


def sample():
    for attempt in range(5):
        before = snapshot()
        free = shutil.disk_usage(BASE).free
        if before == snapshot():
            result = forecast(free, before)
            result['UTC'] = datetime.datetime.now(datetime.timezone.utc).isoformat()
            return result
        time.sleep(0.2)
    return {'status': 'waiting_for_stable_snapshot', 'forecast_all_runs_fit': False,
            'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'all_standard_perambulators_retained': True}


def main(once=False):
    stopped = []
    for sig in [signal.SIGINT, signal.SIGTERM]:
        signal.signal(sig, lambda number, frame: stopped.append(number))
    while not stopped:
        result = sample()
        result['supervisor_PID'] = os.getpid()
        write(WORK / 'batch-storage-ready.json', result)
        write(WORK / 'disk-budget-current.json', result)
        write(WORK / 'storage-supervisor.json',
              {key: value for key, value in result.items() if key not in ('events', 'snapshot')})
        states = result.get('snapshot', [])
        if once or (states and all(state['queue_status'] == 'complete' for state in states)):
            return
        # Sleep in short pieces so an explicit stop is handled promptly.
        for i in range(30):
            if stopped:
                return
            time.sleep(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--once', action='store_true')
    main(parser.parse_args().once)
