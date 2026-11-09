#!/usr/bin/env python3
"""Summarize a completed tuning run only after original Chroma verification."""
import argparse
import csv
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import re
import statistics

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('directory', type=Path)
p.add_argument('--baseline', type=Path, default=Path('results/performance24'))
p.add_argument('--candidates', type=Path, nargs='*', default=[])
a = p.parse_args(); d = a.directory

def rows(path):
    with path.open(newline='') as f:
        return list(csv.DictReader(f))

settings = json.loads((d/'settings.json').read_text())
baseline = json.loads((a.baseline/'summary.json').read_text())
solves = [r for r in rows(d/'test.metal-solves.csv') if r['method'] in ('metal-mg','metal-bicgstab')]
checks = rows(d/'test.chroma-checks.csv')
log = (d/'run.log').read_text()
verify_log = (d/'verify.log').read_text()
if log.count('MG_READY ') != 1 or 'PASS mode=verify' not in verify_log:
    raise RuntimeError('A single hierarchy and successful original Chroma verification are required')
if settings['full_tolerance'] != 1e-12:
    raise RuntimeError('Unexpected full tolerance')
count = settings['rhs']
for method in ['metal-mg','metal-bicgstab']:
    for rhs in range(count):
        measured = [r for r in solves if r['method']==method and int(r['rhs'])==rhs]
        verified = [r for r in checks if r['method']==method and int(r['rhs'])==rhs]
        if len(measured)!=1 or len(verified)!=1:
            raise RuntimeError(f'Missing or duplicate complete solve/check: {method}, RHS {rhs}')
        for r in measured+verified:
            rho = float(r['full_relative'])
            if not math.isfinite(rho) or not 0 <= rho <= 1e-12:
                raise RuntimeError('Full FP64 check failed')
setup = float(re.search(r'MG_READY setup_seconds=([\d.eE+-]+)',log)[1])
selected = re.search(r'SELECTED method=(\S+) pre=(\d+) post=(\d+) weight=([\d.eE+-]+) bottom_cap=(\d+)',log)
if not selected:
    raise RuntimeError('Missing selected cycle')
summary = {}
for method in ['metal-bicgstab','metal-mg']:
    measured = [r for r in solves if r['method']==method]
    verified = [r for r in checks if r['method']==method]
    total = sum(float(r['seconds']) for r in measured)
    hierarchy = setup if method=='metal-mg' else 0
    summary[method] = {
        'rhs_count':count,'solve_seconds':total,'mean_seconds_per_rhs':total/count,
        'hierarchy_setup_seconds':hierarchy,'mean_seconds_per_rhs_with_hierarchy':(total+hierarchy)/count,
        'fine_schur_applications':sum(int(r['fine_schur_applications']) for r in measured),
        'krylov_iterations':sum(int(r['krylov']) for r in measured),
        'full_cpu_fp64_check_seconds':sum(float(r['residual_seconds']) for r in measured),
        'max_native_full_relative':max(float(r['full_relative']) for r in verified),
        'max_relative_solution_difference_from_chroma':max(float(r['relative_solution_difference']) for r in verified),
    }
bi,mg = summary['metal-bicgstab'],summary['metal-mg']
old = baseline['summary']['metal-mg']
comparisons = {
    'fresh_bicgstab_over_tuned_mg_solve_time':bi['solve_seconds']/mg['solve_seconds'],
    'fresh_bicgstab_over_tuned_mg_with_setup':bi['mean_seconds_per_rhs']/mg['mean_seconds_per_rhs_with_hierarchy'],
    'fine_schur_reduction_vs_fresh_bicgstab':1-mg['fine_schur_applications']/bi['fine_schur_applications'],
    'fine_schur_reduction_vs_previous_mg':1-mg['fine_schur_applications']/old['fine_schur_applications'],
    'previous_mg_over_tuned_solve_time':old['mean_seconds_per_rhs']/mg['mean_seconds_per_rhs'],
    'break_even_rhs_for_hierarchy_vs_fresh_bicgstab':math.ceil(setup/(bi['mean_seconds_per_rhs']-mg['mean_seconds_per_rhs'])) if bi['mean_seconds_per_rhs']>mg['mean_seconds_per_rhs'] else None,
}
screens=[]
for candidate in a.candidates:
    params = json.loads((candidate/'settings.json').read_text())
    candidates = rows(candidate/'test.screens.csv')
    outcomes = rows(candidate/'test.metal-solves.csv')
    for r in candidates:
        complete = [v for v in outcomes if v['method']==r['method']]
        screens.append({'directory':str(candidate),'settings':params,'cycle':r,'partial_result':complete})
result = {
    'completed_utc':datetime.now(timezone.utc).isoformat(),'input':baseline['input'],
    'settings':settings,'selected_cycle':{'screen':selected[1],'pre':int(selected[2]),'post':int(selected[3]),'omega':float(selected[4]),'bottom_limit':int(selected[5])},
    'resident_rhs':1,'source_tags':[20261107+i for i in range(count)],
    'screening_source_tag':20261107,'holdout_source_tags':[20261107+i for i in range(1,count)],
    'hierarchies_built_in_final_run':1,
    'timing_scope':baseline['timing_scope'],
    'basis_filter_fine_schur_applications':settings['basis_vectors_per_chirality']*settings['GPU_FP32_filter_sweeps'],
    'previous_baseline_directory':str(a.baseline),'previous_baseline':baseline['summary'],
    'summary':summary,'comparisons':comparisons,'metal_solves':solves,'native_full_checks':checks,
    'screens':screens,
    'screen_note':'Two-correction screens rank candidates; their partial residuals are not accepted full solves. Tuning cost is outside production solve and hierarchy setup timings.',
}
(d/'summary.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({'summary':summary,'comparisons':comparisons},indent=2))
