#!/usr/bin/env python3
"""Audit the new random-source test and compare every complex solution value."""
import array
import csv
import hashlib
import json
import math
import re
import argparse
from datetime import datetime, timezone
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument('--folder', type=Path, default=Path('results/correctness12'))
parser.add_argument('--bridge-prefix', type=Path, default=Path('build/correctness12'))
args = parser.parse_args()
folder = args.folder
prefix = folder/'test'
info = json.loads((folder/'input.json').read_text())
n = math.prod(info['shape'])*12
methods = ['chroma-bicgstab', 'metal-bicgstab', 'metal-mg']
def rows(path):
    with Path(path).open() as f:
        return list(csv.DictReader(f))
checks = rows(str(prefix)+'.chroma-checks.csv')
assert [r['method'] for r in checks] == methods
for r in checks:
    value = float(r['full_relative'])
    assert math.isfinite(value) and 0 <= value <= 1e-12
fields = {}
for method in methods:
    a = array.array('d')
    with Path(str(prefix)+'.'+method+'.bin').open('rb') as f:
        a.fromfile(f, 2*n)
        assert not f.read(1)
    assert all(math.isfinite(x) for x in a)
    fields[method] = a
def compare(a, b):
    norm = math.fsum(x*x for x in b)
    error = math.fsum((x-y)**2 for x, y in zip(a,b))
    maximum = max(math.hypot(a[j]-b[j], a[j+1]-b[j+1]) for j in range(0,len(a),2))
    return {'relative_2norm': math.sqrt(error/norm), 'max_absolute_complex_component': maximum}
pairs = {}
for a,b in [('metal-bicgstab','chroma-bicgstab'),('metal-mg','chroma-bicgstab'),('metal-mg','metal-bicgstab')]:
    pairs[a+'_vs_'+b] = compare(fields[a],fields[b])
    print(a, 'vs', b, pairs[a+'_vs_'+b], flush=True)
native_log = (folder/'chroma-reference.log').read_text()
metal_log = (folder/'metal.log').read_text()
verify_log = (folder/'chroma-verify.log').read_text()
assert 'REFERENCE PASS' in native_log and 'VERIFY PASS' in verify_log and 'METAL PASS' in metal_log
match = lambda pattern, text: float(re.search(pattern, text)[1])
blocks = re.search(r'MG_SETTINGS blocks=(\d+),(\d+),(\d+),(\d+)', metal_log)
fine_blocks = list(map(int, blocks.groups())) if blocks else [3,3,3,4]
summary = {
    'completed_utc': datetime.now(timezone.utc).isoformat(),
    'all_full_fp64_residuals_passed': True,
    'max_native_full_relative': max(float(r['full_relative']) for r in checks),
    'input': info, 'complex_components_per_solution': n,
    'residual_definition': '||b - M x||_2 / ||b||_2 over every site, color and spin of the original full, unpreconditioned Chroma operator',
    'native_full_fp64_checks': checks,
    'full_solution_comparisons': pairs,
    'operator_checks': {
        'native_StoutFermState_vs_exported_links': match(r'link_match=(\S+)', native_log),
        'native_original_full_operator_vs_exported_probe': match(r'UnprecCloverLinOp match=(\S+)', native_log),
        'standalone_original_fp64_vs_native_probe': match(r'original FP64 operator_match=(\S+)', metal_log),
        'metal_vs_cpu_fp32_Schur': match(r'Schur_match=(\S+)', metal_log),
    },
    'native_reference': rows(str(prefix)+'.chroma-reference.csv')[0],
    'metal_solves': rows(str(prefix)+'.metal-solves.csv'),
    'full_fp64_refinement_trace': rows(str(prefix)+'.trace.csv'),
    'mg': {'levels': 3, 'fine_blocks': fine_blocks, 'second_blocks': [2,2,2,2],
           'vectors_per_chirality': [12,12], 'filter_sweeps': 8,
           'filter_weight': .65, 'fine_restart': 16, 'coarse_restart': 8,
           'coarse_limit': 32, 'bottom_limit': 64, 'coarse_target': .05,
           'inner_target': .002, 'full_fp64_target': 1e-12,
           'hierarchy_seconds': match(r'HIERARCHY seconds=(\S+)', metal_log)},
    'precision_and_preconditioning': {
        'reference': 'Original Chroma UnprecCloverLinOp and InvBiCGStab, compiled locally in FP64; target 5e-14, explicitly recomputed full residual checked against 1e-12.',
        'metal_bicgstab': 'FP32 odd asymmetric clover Schur BiCGStab with FP64 full residual refinement and accumulation; zero initial solution.',
        'metal_mg': 'FP32 three-level MG-preconditioned FGMRES on the same odd asymmetric Schur operator, with the same FP64 refinement; zero initial solution.',
        'schur': 'A_oo - H_oe A_ee^-1 H_eo, with -1/2 included in H.',
        'note': 'Residual tolerance and relative solution difference are distinct quantities. The reference uses the full unpreconditioned system; iteration counts are not an apples-to-apples speed comparison.',
    },
    'native_build': {
        'chroma_config': '/Users/admin/Downloads/macos-chroma/install/chroma/bin/chroma-config',
        'matching_source': '/Users/admin/Downloads/qcd-macos-share-20261005/chroma.devel',
        'BASE_PRECISION': 64,
        'note': 'Original native parameter parser, stout state parameters, stout kernels, UnprecCloverLinOp and InvBiCGStab compiled locally in FP64. External Chroma/SuperbBLAS trees unchanged.',
    },
}
if (folder/'workflow.json').exists():
    summary['workflow'] = json.loads((folder/'workflow.json').read_text())
code = sorted(set(Path('.').glob('*.hpp')) | set(Path('.').glob('*.metal')) |
              {Path(p) for p in ['correctness12_prepare.py','correctness12_report.py','correctness_chroma.cc',
                                'correctness_compare.mm','production_export.cc','Makefile.production','README.md']})
files = [Path(info['gauge']), Path(str(args.bridge_prefix)+'.bridge'), Path(str(args.bridge_prefix)+'.rhs.bin'),
         Path('build/correctness_chroma'),Path('build/correctness_compare')] + code
files += [p for p in folder.iterdir() if p.is_file() and p.name not in ('summary.json','files.sha256')]
summary['file_hashes'] = []
for path in files:
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for data in iter(lambda: stream.read(8*1024*1024), b''):
            h.update(data)
    summary['file_hashes'].append({'path':str(path.resolve()), 'bytes':path.stat().st_size, 'sha256':h.hexdigest()})
(folder/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
(folder/'files.sha256').write_text(''.join(r['sha256']+'  '+r['path']+'\n' for r in summary['file_hashes']))
print('AUDIT PASS: all three independent full FP64 residuals <= 1e-12; compared all', n, 'complex components', flush=True)
