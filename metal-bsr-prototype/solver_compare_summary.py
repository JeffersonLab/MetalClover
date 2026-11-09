#!/usr/bin/env python3
"""Audit paired RHS timings and compare all six complete solver paths."""
import argparse
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('csv',nargs='?',type=Path,default=Path('results/solver-compare16.csv'))
ap.add_argument('--output',type=Path,default=Path('SOLVER_COMPARISON_RESULTS.md'))
args=ap.parse_args()
with args.csv.open() as f: rows=list(csv.DictReader(f))
with Path(str(args.csv)+'.trace.csv').open() as f: traces=list(csv.DictReader(f))
groups=defaultdict(list)
for r in rows:
    assert int(r['extent'])==16 and float(r['csw'])==1
    assert int(r['tile_rhs'])==96 and int(r['rhs'])==96*int(r['tiles'])
    assert int(r['bi_chunk'])==1 and float(r['inner_tolerance'])==.002
    assert all(math.isfinite(float(r[k])) for k in ('total_ms','max_true_residual','solution_error','worst_coarse_true_residual'))
    assert float(r['max_true_residual'])<=float(r['tolerance'])==1e-10
    assert float(r['solution_error'])<=1e-7
    if not r['backend'].endswith('sparse_mg'):
        assert float(r['hierarchy_setup_ms'])==0 and int(r['coarse_spmv_calls'])==0
    groups[float(r['mass']),int(r['rhs']),r['backend']].append(r)
cases=sorted({k[:2] for k in groups})
methods=('bicgstab','smooth','sparse_mg')
backends=tuple(p+'_'+m for p in ('cpu','metal') for m in methods)
assert cases==[(.1,96),(.1,192),(.5,96),(.5,192)]
for k in cases:
    for b in backends: assert len(groups[k+(b,)])==3,(k,b)
assert len(rows)==72
last=defaultdict(dict)
for t in traces:
    key=(float(t['mass']),t['backend'],int(t['trial']))
    step,offset=int(t['outer_step']),int(t['offset'])
    assert math.isfinite(float(t['max_true_residual']))
    if offset not in last[key] or step>int(last[key][offset]['outer_step']):last[key][offset]=t
for r in rows:
    final={offset:t for offset,t in last[float(r['mass']),r['backend'],int(r['trial'])].items() if offset<int(r['rhs'])}
    assert len(final)==int(r['tiles'])
    assert all(float(t['max_true_residual'])<=1e-10 for t in final.values())
    assert math.isclose(max(float(t['max_true_residual']) for t in final.values()),float(r['max_true_residual']),rel_tol=1e-10)
physical=sum(len(v) for v in last.values())
assert physical==72

def values(k,b,field):return [float(r[field]) for r in groups[k+(b,)]]
def med(k,b,field='total_ms'):return st.median(values(k,b,field))
def span(k,b,field):
    v=sorted({int(r[field]) for r in groups[k+(b,)]})
    return str(v[0]) if len(v)==1 else f'{v[0]}–{v[-1]}'
def winner(k,p):return min(methods,key=lambda m:med(k,p+'_'+m))
labels={'bicgstab':'BiCGStab','smooth':'smoother FGMRES','sparse_mg':'sparse MG'}
lines=['# BiCGStab, smoother FGMRES, and sparse multigrid on 16⁴','',
       'Fresh comparison on the Apple M1 Max (24 GPU cores, 32 GB, macOS 26.7.1). '
       'All methods use the same 16⁴ synthetic clover-Wilson operator, c_sw=1, '
       'gauge strength 0.25, QDP spin conventions, antiperiodic time, and positive '
       'bare masses 0.1/0.5. All RHS, odd Schur reduction/reconstruction, and '
       'stopping targets agree across methods. Inner work is FP32; accumulation '
       'and every component of the full reconstructed residual are CPU FP64. '
       'Every RHS must meet a true relative residual of 1e-10.','',
       '## Complete solve times','',
       'Medians of three trials from zero, after one complete resident-tile '
       'warmup per method. **96 RHS are resident. 192 RHS use two sequential '
       '96-RHS tiles.** Within a trial the first measured tile supplies the '
       '96-RHS time, and both actual measured tiles supply the 192-RHS total. '
       'Thus the two sizes are paired prefix measurements; 192 is not an '
       'extrapolation and is not 192 simultaneously resident RHS.','',
       '| Mass | RHS | CPU BiCG s | CPU smoother FG s | CPU sparse MG s | Metal BiCG s | Metal smoother FG s | Metal sparse MG s |',
       '|---:|---:|---:|---:|---:|---:|---:|---:|']
for k in cases:
    lines.append(f'| {k[0]:g} | {k[1]} | '+' | '.join(f'{med(k,b)/1000:.3f}' for b in backends)+' |')
lines+=['','Timings include zero/reset work, FP32 preparation/reconstruction, '
        'all fine and coarse products, vector operations/reductions, small '
        'CPU QR, GPU submissions/waits, FP64 solution updates and full FP64 '
        'residual checks. Setup and manufactured RHS construction are outside '
        'solve time and reported separately. Known-solution error is checked '
        'after acceptance outside timing. CPU uses eight workers. Methods '
        'run sequentially, rotating order across mass/processor cases, with '
        'processor order reversed between masses. Ordinary local system '
        'activity was present; thermal and scheduling state were not isolated.','',
        '| Mass | RHS | Fastest CPU median | Fastest Metal median | CPU BiCG / MG | CPU smoother / MG | Metal BiCG / MG | Metal smoother / MG |',
        '|---:|---:|---|---|---:|---:|---:|---:|']
for k in cases:
    gains=[med(k,p+'_'+m)/med(k,p+'_sparse_mg') for p in ('cpu','metal') for m in ('bicgstab','smooth')]
    lines.append(f'| {k[0]:g} | {k[1]} | {labels[winner(k,"cpu")]} | {labels[winner(k,"metal")]} | '+' | '.join(f'{v:.2f}×' for v in gains)+' |')
lines+=['','A baseline/MG ratio above one favors sparse MG; below one favors '
        'the baseline. Fastest refers to the observed median. Small '
        'differences within the timing spread do not establish a robust ranking.','',
        '| Mass | RHS | BiCG CPU / Metal | Smoother FG CPU / Metal | Sparse MG CPU / Metal |',
        '|---:|---:|---:|---:|---:|']
for k in cases:
    lines.append(f'| {k[0]:g} | {k[1]} | '+' | '.join(f'{med(k,"cpu_"+m)/med(k,"metal_"+m):.2f}×' for m in methods)+' |')
lines+=['','## Convergence and actual operator work','',
        'BiCGStab and FGMRES iterations represent different work. Fine S '
        'counts include smoothing, Arnoldi products, true FP32 refreshes, '
        'and every dispatched BiCGStab product. Each call applies to one '
        'resident 96-RHS tile; 192-RHS rows sum both tiles. GPU BiCGStab '
        'checks convergence every iteration (chunk=1), minimizing iteration '
        'padding; its second half-step is still dispatched when RHS are '
        'already converged. Krylov counts sum the maximum per-RHS iteration '
        'count over outer corrections; ranges cover resident tiles/trials. '
        'All methods retain full FP64 outer checks.','',
        '| Mass | RHS | Method | Outer corrections/tile | Krylov steps/tile | Total fine S | Coarse SpMV | Worst full FP64 residual |',
        '|---:|---:|---|---:|---:|---:|---:|---:|']
for k in cases:
    for b in backends:
        omin,omax=min(values(k,b,'outer_min')),max(values(k,b,'outer_max'))
        kmin,kmax=min(values(k,b,'krylov_min')),max(values(k,b,'krylov_max'))
        outer=f'{omin:g}' if omin==omax else f'{omin:g}–{omax:g}'
        kry=f'{kmin:g}' if kmin==kmax else f'{kmin:g}–{kmax:g}'
        lines.append(f'| {k[0]:g} | {k[1]} | {b} | {outer} | {kry} | {span(k,b,"fine_schur_calls")} | {span(k,b,"coarse_spmv_calls")} | {max(values(k,b,"max_true_residual")):.3g} |')
lines+=['','Smoother FGMRES uses two clover-block Jacobi steps at weight '
        '1.0, with restart 8. Sparse MG uses one pre/post step at the same '
        'weight and restart, 12 basis vectors per chirality, 4⁴ aggregates, '
        'and eight FP64 filtering sweeps at setup weight 0.65. Its 6,144-'
        'unknown coarse operator has 8,448 dense 24×24 sparse blocks. '
        'Coarse BiCGStab has local block preconditioning, target 0.1 in '
        'the preconditioned residual and a four-iteration cap. Flexible '
        'GMRES permits this inexact, variable correction. The actual '
        'unpreconditioned FP32 coarse residual is recomputed and reported '
        'in raw data; it may exceed 0.1. It never replaces the full FP64 '
        'stopping test. GPU coarse counts include masked/padded dispatches.','',
        'Full CPU FP64 history for the second 96-RHS tile (global offset 96), '
        'mass 0.1, trial 0:','',
        '| Outer correction | CPU BiCG | CPU smoother FG | CPU MG | Metal BiCG | Metal smoother FG | Metal MG |',
        '|---:|---:|---:|---:|---:|---:|---:|']
history=defaultdict(dict)
for t in traces:
    if float(t['mass'])==.1 and t['trial']=='0' and t['offset']=='96':history[int(t['outer_step'])][t['backend']]=float(t['max_true_residual'])
for step,v in sorted(history.items()):
    lines.append(f'| {step} | '+' | '.join(f'{v[b]:.3g}' if b in v else 'accepted earlier' for b in backends)+' |')
lines+=['','## Setup and an unlimited supply of RHS','',
        'Each mass-specific hierarchy is built once and reused across all '
        'methods, processors and trials. Baselines need no hierarchy; their '
        'hierarchy setup field is zero even though the shared experiment '
        'also prepares MG. Fine-operator/problem generation and FP64 outer '
        'storage are common. Each solver allocates its own inner workspace '
        'sequentially. Allocation timing does not include all first-use '
        'page/pipeline costs, which warmup covers.','',
        '| Mass | RHS | MG hierarchy s | CPU best baseline | CPU hierarchy-only break-even batches | Metal best baseline | Metal hierarchy-only break-even batches |',
        '|---:|---:|---:|---|---:|---|---:|']
for k in cases:
    setup=med(k,'cpu_sparse_mg','hierarchy_setup_ms')
    cells=[]
    for p in ('cpu','metal'):
        b=min(('bicgstab','smooth'),key=lambda m:med(k,p+'_'+m));saved=med(k,p+'_'+b)-med(k,p+'_sparse_mg')
        cells.extend([labels[b],f'{setup/saved:.1f}' if saved>0 else 'No warm gain'])
    lines.append(f'| {k[0]:g} | {k[1]} | {setup/1000:.3f} | '+' | '.join(cells)+' |')
lines+=['','The estimate divides hierarchy setup by measured warm batch time '
        'saved relative to the faster baseline, excluding method-specific '
        'workspace/compilation differences. It is sensitive to small timing '
        'differences. More RHS can amortize setup when MG has a warm-time '
        'gain; it cannot remove a per-batch slowdown.','',
        '| Mass | Method | Inner workspace setup s | Common outer setup s | Common problem setup s | Peak process RSS GiB |',
        '|---:|---|---:|---:|---:|---:|']
for mass in (.1,.5):
    k=(mass,96)
    for b in backends:
        lines.append(f'| {mass:g} | {b} | {med(k,b,"inner_workspace_setup_ms")/1000:.3f} | {med(k,b,"outer_workspace_setup_ms")/1000:.3f} | {med(k,b,"problem_setup_ms")/1000:.3f} | {max(values(k,b,"peak_rss_mib"))/1024:.2f} |')
lines+=['','Process RSS is the cumulative high-water mark within the run, '
        'rather than an isolated per-method memory measurement. The sparse '
        'FP32 coefficient matrix occupies 37.125 MiB, with a separate FP64 '
        'copy and basis/local inverses. Large fine Krylov vectors dominate '
        'memory, hence sequential method allocation and 96-wide tiling.','',
        '## Timing spread and FP64 cost','',
        '| Mass | RHS | Method | Min–max solve s | Median full FP64 residual s | Median FP32 inner s |',
        '|---:|---:|---|---:|---:|---:|']
for k in cases:
    for b in backends:
        v=values(k,b,'total_ms')
        lines.append(f'| {k[0]:g} | {k[1]} | {b} | {min(v)/1000:.3f}–{max(v)/1000:.3f} | {med(k,b,"residual_ms")/1000:.3f} | {med(k,b,"inner_ms")/1000:.3f} |')
lines+=['','## Validation and reproduction','',
        f'All {len(rows)} paired logical batch records passed every-RHS full '
        f'CPU FP64 acceptance. They derive from {physical} physically '
        f'measured 96-RHS tile solves ({physical*96:,} RHS executions), '
        'with shared first tiles between the reported logical sizes. '
        f'Worst final relative residual: {max(float(r["max_true_residual"]) for r in rows):.3g}; '
        f'maximum known-solution error: {max(float(r["solution_error"]) for r in rows):.3g}. '
        'Every CSV acceptance value was audited against the final per-tile '
        'FP64 trace values. The 48-record small comparison self-test passed '
        'all six paths, two volumes and masses, and zero RHS. The 40-solve '
        'established clover regression passed after ensuring GPU BiCGStab '
        'allocates a full direction workspace when optional multiply '
        'inputs have been released. Sparse Galerkin validation from the '
        'previous study remains applicable: no hierarchy or fine/coarse '
        'operator algorithm changed.','',
        '```sh','make build/metal_solver_compare',
        'VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare --self-test',
        'VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \\',
        '  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass .1,.5 --csw 1 \\',
        '  --samples 3 --threads 8 --bi-chunk 1 --csv results/solver-compare16.csv',
        'python3 solver_compare_summary.py','```','',
        'See [SOLVER_COMPARISON.md](SOLVER_COMPARISON.md) for algorithms and '
        'the paired timing protocol. Results describe synthetic positive-'
        'mass fields on this M1 Max, not a production ensemble or critical-'
        'mass scaling study. Chroma and SuperbBLAS sources are unchanged.','',
        f'Raw [measurements]({args.csv.as_posix()}), '
        f'[full FP64 traces]({args.csv.as_posix()}.trace.csv), '
        '[run log](results/solver-compare16.log).','']
args.output.write_text('\n'.join(lines))
print(f'Wrote {args.output}: {len(rows)} records / {physical} physical tile solves')
