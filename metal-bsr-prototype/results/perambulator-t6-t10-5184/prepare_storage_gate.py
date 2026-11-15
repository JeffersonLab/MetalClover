"""Forecast all remaining outputs before enabling the detached batch queue.

This preparation script requires the observed time-four solve/time-five-wait
state. It never deletes files or starts an inversion.
"""
from pathlib import Path
import datetime
import json
import shutil

REPORT = Path(__file__).resolve().parent
BASE = Path('/Users/admin/Downloads/MetalClover-results')
GATE = Path('/private/tmp/chroma-times6-10-production/batch-storage-ready.json')
GIB = 2**30
POINT = 6879707264
RAW_MATRIX = 16 * 5184**2 * 16
STANDARD_BUDGET = RAW_MATRIX + 64 * 2**20
PER_JOB_EVIDENCE_RESERVE = GIB // 4
OTHER_RESERVE = GIB


def forecast(free_bytes, point4_allocated):
    # Allow for Git records and unrelated modest disk activity up front.
    free = free_bytes - OTHER_RESERVE
    events = [{'event': 'global 1 GiB reserve', 'free_GiB': free / GIB}]
    starts, projection_starts = [], []
    kept_points = {4: point4_allocated}
    for t in range(4, 11):
        if t == 6:
            free += sum(kept_points.values())
            kept_points.clear()
            events.append({'event': 'retire accepted points 4/5 after full evidence',
                           'free_GiB': free / GIB})
        if t >= 6:
            starts.append(free / GIB)
        events.append({'event': f'time-{t} start', 'free_GiB': free / GIB})
        if t == 4:
            free -= max(0, POINT - point4_allocated)
        else:
            free -= POINT
            kept_points[t] = POINT
        free -= PER_JOB_EVIDENCE_RESERVE
        projection_starts.append(free / GIB)
        events.append({'event': f'time-{t} before standard projection',
                       'free_GiB': free / GIB})
        free -= STANDARD_BUDGET
        events.append({'event': f'time-{t} peak output allocation',
                       'free_GiB': free / GIB})
        if t >= 6:
            free += kept_points.pop(t)
            events.append({'event': f'time-{t} accepted point retired',
                           'free_GiB': free / GIB})
    minimum = min(event['free_GiB'] for event in events)
    fits = (minimum >= 4.0 and min(starts) >= 18.0
            and min(projection_starts) >= (1.1 * RAW_MATRIX / GIB + 4.0))
    return {'status': 'pass' if fits else 'insufficient_space',
            'forecast_all_runs_fit': fits,
            'all_standard_perambulators_retained': True,
            'minimum_forecast_free_GiB': minimum,
            'minimum_new_job_start_free_GiB': min(starts),
            'minimum_before_projection_GiB': min(projection_starts),
            'native_projection_required_free_GiB': 1.1 * RAW_MATRIX / GIB + 4.0,
            'final_forecast_free_GiB': free / GIB,
            'standard_output_budget_bytes_each': STANDARD_BUDGET,
            'per_job_log_and_evidence_reserve_bytes': PER_JOB_EVIDENCE_RESERVE,
            'additional_global_reserve_bytes': OTHER_RESERVE,
            'events': events}


def main():
    assert not GATE.exists(), 'Never overwrite an active storage gate'
    point4 = BASE / 'cfg1000-t4-unsmeared.npy'
    pipeline = json.loads((BASE / 'cfg1000-t4-unsmeared.pipeline.json').read_text())
    assert pipeline['status'] == 'running' and pipeline['stage'] == 'solve'
    assert point4.stat().st_size == POINT and not point4.is_symlink()
    for t in range(4, 11):
        standard = BASE / f'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb1000'
        assert not standard.exists()
        if t >= 5:
            assert not (BASE / f'cfg1000-t{t}-unsmeared.npy').exists()
            assert not (BASE / f'cfg1000-t{t}-unsmeared.pipeline.json').exists()
    preservation = json.loads((REPORT / 'storage-preservation.json').read_text())
    assert preservation['status'] == 'complete'
    binaries = json.loads((REPORT / 'binary-preservation/storage-preservation.json').read_text())
    assert binaries['status'] == 'complete'
    free = shutil.disk_usage(BASE).free
    record = forecast(free, point4.stat().st_blocks * 512)
    record.update(UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  actual_free_disk_GiB=free / GIB,
                  initial_state='time four solving with its full point file allocated; time five waiting',
                  storage_preservation_status='complete',
                  assumes_no_unbudgeted_large_disk_writes=True)
    (REPORT / 'disk-budget.json').write_text(json.dumps(record, indent=2) + '\n')
    assert record['forecast_all_runs_fit'], record
    temporary = GATE.with_suffix('.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    temporary.replace(GATE)
    print(json.dumps({key: value for key, value in record.items() if key != 'events'}, indent=2))


if __name__ == '__main__':
    main()
