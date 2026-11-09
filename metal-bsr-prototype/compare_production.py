#!/usr/bin/env python3
"""Audit the four v0 solves and compare the stored perambulator elements."""
import argparse
import csv
import json
import math
import re
from pathlib import Path

parser = argparse.ArgumentParser()
parser.add_argument('prefix', type=Path)
parser.add_argument('--reference', type=Path, default=Path('results/anisotropic24-peram-00.csv'))
parser.add_argument('--source-time', type=int, default=176)
parser.add_argument('--first-time', type=int, default=176)
parser.add_argument('--times', type=int, default=48)
parser.add_argument('--chroma-residuals', type=Path,
                    help='independent native Chroma FP64 residual CSV')
parser.add_argument('--solve-log', type=Path, help='include measured solve/setup times')
parser.add_argument('--action', type=Path, help='include the action and input provenance JSON')
args = parser.parse_args()
with args.reference.open() as stream:
    reference = list(csv.DictReader(stream))
with Path(str(args.prefix)+'.csv').open() as stream:
    computed = list(csv.DictReader(stream))
with Path(str(args.prefix)+'.trace.csv').open() as stream:
    traces = list(csv.DictReader(stream))

expected = {(s, k, t) for s in range(4) for k in range(4)
            for t in range(args.first_time, args.first_time+args.times)}
ref = {}
for row in reference:
    if int(row['t_source']) != args.source_time:
        continue
    key = (int(row['spin_source']), int(row['spin_sink']), int(row['t_sink']))
    if key in ref:
        raise ValueError('duplicate reference element')
    ref[key] = complex(float(row['real']), float(row['imag']))
got = {}
for row in computed:
    key = (int(row['spin_source']), int(row['spin_sink']), int(row['time']))
    if key in got:
        raise ValueError('duplicate computed element')
    got[key] = complex(float(row['real']), float(row['imag']))
if ref.keys() != expected or got.keys() != expected:
    raise ValueError('perambulator element coverage differs from the requested four-spin window')
if not all(math.isfinite(z.real) and math.isfinite(z.imag)
           for values in (ref, got) for z in values.values()):
    raise ValueError('nonfinite perambulator element')

final = {}
for row in traces:
    spin = int(row['spin_source'])
    value = float(row['full_relative'])
    if not math.isfinite(value):
        raise ValueError('nonfinite full residual trace')
    final[spin] = value
if set(final) != set(range(4)) or max(final.values()) > 1e-10:
    raise ValueError('four accepted full FP64 residuals are required')

native = None
if args.chroma_residuals:
    with args.chroma_residuals.open() as stream:
        rows = list(csv.DictReader(stream))
    native = {int(r['spin_source']): float(r['full_relative']) for r in rows}
    if (len(rows) != 4 or set(native) != set(range(4)) or
            not all(math.isfinite(v) and 0 <= v <= 1e-10 for v in native.values())):
        raise ValueError('four accepted independent Chroma FP64 residuals are required')

def metric(keys):
    norm = sum(abs(ref[k])**2 for k in keys)
    error = sum(abs(got[k]-ref[k])**2 for k in keys)
    worst = max(keys, key=lambda k: abs(got[k]-ref[k]))
    return {
        'relative_2norm': math.sqrt(error/norm),
        'max_absolute': abs(got[worst]-ref[worst]),
        'max_absolute_at': {'spin_source': worst[0], 'spin_sink': worst[1], 'time': worst[2]},
    }

summary = {
    'elements': len(expected),
    'source_time': args.source_time,
    'sink_times': [args.first_time, args.first_time+args.times-1],
    'final_full_residual_by_source_spin': final,
    'all_elements': metric(expected),
    'by_source_spin': {s: metric([k for k in expected if k[0] == s]) for s in range(4)},
    'by_sink_time': {t: metric([k for k in expected if k[2] == t])
                     for t in range(args.first_time, args.first_time+args.times)},
}
example = (0, 0, args.source_time)
if example in expected:
    summary['same_time_spin00_example'] = {
        'computed': [got[example].real, got[example].imag],
        'reference': [ref[example].real, ref[example].imag],
        'absolute_difference': abs(got[example]-ref[example]),
    }
if native is not None:
    summary['independent_chroma_full_residual_by_source_spin'] = native
if args.action:
    summary['action'] = json.loads(args.action.read_text())
if args.solve_log:
    log = args.solve_log.read_text()
    measured = {}
    for s, seconds, residual_s, k, residual in re.findall(
            r'SOLVED spin=(\d+) seconds=(\S+) residual_s=(\S+) K=(\d+) full_relative=(\S+)', log):
        measured[int(s)] = {'seconds': float(seconds), 'fp64_residual_seconds': float(residual_s),
                            'krylov_iterations': int(k), 'final_full_relative': float(residual)}
    if set(measured) != set(range(4)):
        raise ValueError('complete solve timing records are required')
    hierarchy = re.search(r'HIERARCHY built seconds=(\S+) coarse1=(\d+) coarse2=(\d+)', log)
    if not hierarchy:
        raise ValueError('missing hierarchy record')
    summary['timing'] = {
        'per_source_spin': measured,
        'hierarchy_seconds': float(hierarchy[1]),
        'coarse_unknowns': [int(hierarchy[2]), int(hierarchy[3])],
        'four_solve_seconds': sum(m['seconds'] for m in measured.values()),
        'fp64_residual_seconds': sum(m['fp64_residual_seconds'] for m in measured.values()),
        'initial_guess_spins': [int(s) for s in re.findall(r'INITIAL guess spin=(\d+)', log)],
        'note': 'Solve times include mixed precision corrections and CPU FP64 residuals; '
                'LIME import, stout smearing, setup, checkpoint I/O and native verification are separate. '
                'An initial guess is used only for the spins explicitly listed.',
    }
with Path(str(args.prefix)+'.comparison.csv').open('w', newline='') as stream:
    writer = csv.writer(stream)
    writer.writerow(['spin_source', 'spin_sink', 'time', 'computed_real',
                     'computed_imag', 'reference_real', 'reference_imag',
                     'difference_real', 'difference_imag', 'absolute_difference'])
    for s, k, t in sorted(expected):
        a, b = got[s, k, t], ref[s, k, t]
        d = a-b
        writer.writerow([s, k, t, a.real, a.imag, b.real, b.imag,
                         d.real, d.imag, abs(d)])
Path(str(args.prefix)+'.comparison.json').write_text(json.dumps(summary, indent=2)+'\n')
print(json.dumps({k: v for k, v in summary.items() if k != 'by_sink_time'}, indent=2))
