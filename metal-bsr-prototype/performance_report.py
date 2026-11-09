#!/usr/bin/env python3
"""Summarize completed real-field timings, requiring all native FP64 checks."""
import argparse
import csv
from datetime import datetime, timezone
import json
from pathlib import Path
import re

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('directory', type=Path)
p.add_argument('--rhs', type=int, required=True)
a = p.parse_args(); d = a.directory
def rows(name):
    with (d/name).open() as f:
        return list(csv.DictReader(f))
native = rows('test.chroma-solves.csv')
metal = rows('test.metal-solves.csv')
checks = rows('test.chroma-checks.csv')
methods = ['chroma-bicgstab', 'metal-bicgstab', 'metal-mg']
expected = {(str(rhs), method) for rhs in range(a.rhs) for method in methods}
assert {(r['rhs'], r['method']) for r in checks} == expected
assert len(checks) == 3*a.rhs
assert all(float(r['full_relative']) <= 1e-12 for r in checks)
assert len(native) == a.rhs and len(metal) == 2*a.rhs
assert {r['rhs'] for r in native} == {str(i) for i in range(a.rhs)}
assert {(r['rhs'],r['method']) for r in metal} == expected-{(str(i),'chroma-bicgstab') for i in range(a.rhs)}
log = (d/'metal.log').read_text()
blocks = list(map(int,re.search(r'MG_SETTINGS blocks=([\d,]+)',log)[1].split(',')))
setup = float(re.search(r'MG_READY setup_seconds=([\d.e+-]+)', log)[1])
allocated = float(re.search(r'allocated_GiB=([\d.e+-]+)', log)[1])
cpu_log = (d/'chroma.log').read_text()
cpu_setup = float(re.search(r'setup_seconds=([\d.e+-]+)', cpu_log)[1])
assert 'PASS mode=verify' in (d/'verify.log').read_text()
state_log = d/'stout-state-check.log'
state_match = None
if state_log.exists():
    state_text = state_log.read_text()
    assert 'STATE PASS' in state_text
    state_match = float(re.search(r'link_match=([\d.e+-]+)',state_text)[1])
    assert state_match<=1e-14
summary = {}
for method in methods:
    r = native if method == methods[0] else [r for r in metal if r['method']==method]
    total = sum(float(x['seconds']) for x in r)
    method_setup = setup if method == 'metal-mg' else 0
    summary[method] = {
        'rhs_count': a.rhs,
        'solve_seconds': total,
        'mean_seconds_per_rhs': total/a.rhs,
        'hierarchy_setup_seconds': method_setup,
        'solve_plus_hierarchy_seconds': total+method_setup,
        'mean_seconds_per_rhs_with_hierarchy': (total+method_setup)/a.rhs,
        'max_native_full_relative': max(float(x['full_relative']) for x in checks if x['method']==method),
        'max_relative_solution_difference_from_chroma': max(float(x['relative_solution_difference']) for x in checks if x['method']==method),
    }
cpu_time = summary['chroma-bicgstab']['solve_seconds']
for method in methods[1:]:
    summary[method]['cpu_over_metal_solve_time'] = cpu_time/summary[method]['solve_seconds']
    summary[method]['cpu_over_metal_solve_plus_hierarchy_time'] = cpu_time/summary[method]['solve_plus_hierarchy_seconds']
result = {
    'completed_utc': datetime.now(timezone.utc).isoformat(),
    'input': json.loads((d/'input.json').read_text()),
    'resident_rhs': 1,
    'source_tags': list(range(20261107,20261107+a.rhs)),
    'source_initial_guess': 'zero for every solver and RHS',
    'full_residual_definition': '||b-Mx||_2/||b||_2 over all sites, spins, and colors, using the original Chroma FP64 operator',
    'native_cpu': 'Original Chroma InvBiCGStabReliable, original FP64 EvenOddPrecCloverLinOp, original FP32 EvenOddPrecDumbCloverFLinOp; Delta=0.1, Schur target=5e-14, final full target=1e-12',
    'metal': 'FP32 even/odd BiCGStab or MG-preconditioned FGMRES, CPU FP64 full defect correction with inner target 0.002 and final target 1e-12',
    'mg': {'levels': 3, 'fine_aggregation': blocks, 'second_aggregation': [2,2,2,2], 'vectors_per_chirality_each_level': 12, 'filter_sweeps': 8, 'filter_tile': 2, 'fine_restart': 8, 'first_coarse_restart': 8, 'first_coarse_iteration_cap': 16, 'bottom_iteration_cap': 16, 'coarse_relative_tolerance': 0.1, 'hierarchies_built': 1},
    'timing_scope': 'Solve includes RHS reduction, reconstruction, vector updates, and every full CPU FP64 residual check. Source generation, operator setup, import, Metal shader compilation, saved-field I/O, and subsequent independent Chroma verification are outside solve timings. Hierarchy setup is reported separately and amortized over the tested RHS.',
    'cpu_operator_import_setup_seconds': cpu_setup,
    'metal_allocated_gib_at_mg_ready': allocated,
    'native_stout_state_link_match': state_match,
    'summary': summary,
    'native_solves': native,
    'metal_solves': metal,
    'native_full_checks': checks,
}
(d/'summary.json').write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(summary, indent=2))
