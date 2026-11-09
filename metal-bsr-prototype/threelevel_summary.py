#!/usr/bin/env python3
"""Audit complete physical RHS sweeps, then report measured paired totals."""
import argparse
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('csv',nargs='?',type=Path,default=Path('results/threelevel32-compare.csv'))
ap.add_argument('--output',type=Path,default=Path('THREE_LEVEL_RESULTS.md'))
ap.add_argument('--resumed-after-reboot',action='store_true',
                help='Record the reboot and resumed warmup in the timing protocol.')
args=ap.parse_args()
with args.csv.open() as f:rows=list(csv.DictReader(f))
with Path(str(args.csv)+'.trace.csv').open() as f:traces=list(csv.DictReader(f))
assert rows
masses={float(r['mass']) for r in rows};tiles={int(r['tile_rhs']) for r in rows}
bottoms={int(r['bottom_limit']) for r in rows};samples={int(r['trial']) for r in rows}
assert len(masses)==len(tiles)==len(bottoms)==1
mass=next(iter(masses));tile=next(iter(tiles));bottom=next(iter(bottoms))
assert tile in (4,8) and mass<0 and samples==set(range(len(samples)))
methods=('bicgstab','smooth','mg3');processors=('cpu','metal')
groups=defaultdict(list);last={}
for t in traces:
    assert int(t['extent'])==32 and float(t['mass'])==mass and int(t['tile_rhs'])==tile
    key=(t['backend'],t['method'],int(t['trial']),int(t['offset']))
    if key not in last or int(t['outer'])>int(last[key]['outer']):last[key]=t
keys=set()
for r in rows:
    assert int(r['extent'])==32 and float(r['csw'])==1
    assert int(r['workers'])==8 and int(r['vectors1'])==int(r['vectors2'])==12
    assert int(r['block1'])==4 and int(r['restart'])==8 and int(r['coarse_limit'])==8
    assert float(r['tolerance'])==1e-10 and float(r['inner_tolerance'])==.002
    assert r['method'] in methods and r['backend'] in processors
    assert all(math.isfinite(float(r[k])) for k in ('total_ms','inner_ms','residual_ms','max_true_residual','solution_error','worst_bottom_true_residual'))
    assert float(r['max_true_residual'])<=1e-10 and float(r['solution_error'])<=1e-7
    key=(r['backend'],r['method'],int(r['trial']),int(r['offset']))
    assert key not in keys;keys.add(key)
    assert key in last and math.isclose(float(last[key]['max_true_residual']),float(r['max_true_residual']),rel_tol=1e-11)
    if r['method']!='mg3':assert int(r['level1_products'])==int(r['level2_products'])==0 and float(r['hierarchy_ms'])==0
    groups[r['backend'],r['method'],int(r['trial'])].append(r)
for b in processors:
    for m in methods:
        for trial in samples:assert sorted(int(r['offset']) for r in groups[b,m,trial])==list(range(0,192,tile))
assert len(rows)==6*len(samples)*192//tile
def selected(b,m,trial,rhs):return [r for r in groups[b,m,trial] if int(r['offset'])<rhs]
def total(b,m,rhs,field='total_ms'):
    return st.median(sum(float(r[field]) for r in selected(b,m,t,rhs)) for t in samples)
def worst(b,m,rhs,field):return max(float(r[field]) for t in samples for r in selected(b,m,t,rhs))
labels={'bicgstab':'BiCGStab','smooth':'Smoother FGMRES','mg3':'Three-level MG'}
lines=['# Three-level clover-Wilson comparison at 32⁴','',
       f'Apple M1 Max (24 GPU cores, 32 GB), eight CPU workers. Bare mass **m₀={mass:g}**, '
       '**c_sw=1**, the same synthetic gauge generator (strength 0.25), QDP gamma '
       'conventions, and antiperiodic time. All inner Krylov vectors and '
       'GPU operations are FP32. Every accepted full solution passes a CPU FP64 '
       'residual check of the original reconstructed clover-Wilson equation for '
       'every RHS, with relative 2-norm target 1e-10. CPU fine BiCGStab uses '
       'FP32 block sums and a final FP64 reduction of the partials, rounded '
       'back to FP32 scalar values; Metal reductions stay FP32.','',
       '## Actual total solve time','',
       f'{len(samples)} complete paired sweep(s) per solver/processor. '
       f'**{tile} RHS are resident**. The 96-RHS result sums the first {96//tile} '
       f'measured tiles; the 192-RHS result sums all {192//tile} distinct tiles. '
       'Every tile starts from zero. These are measured totals, not projected '
       'times. One resident-tile warmup precedes each method.','',
       '| Solver | CPU 96 RHS s | Metal 96 RHS s | CPU 192 RHS s | Metal 192 RHS s |',
       '|---|---:|---:|---:|---:|']
for m in methods:lines.append('| '+labels[m]+' | '+' | '.join(f'{total(b,m,n)/1000:.2f}' for n in (96,192) for b in processors)+' |')
lines+=['','Complete solve time includes initialization, elimination and '
        'reconstruction, fine/coarse products and reductions, CPU small QR, '
        'GPU submissions/waits, FP64 updates, full CPU FP64 residual checks, '
        'and per-outer progress logging. Setup, manufactured RHS generation, '
        'warmup, and known-solution error checking are outside solve time. '
        'Methods run sequentially, Metal first, under ordinary desktop load.','']
if args.resumed_after_reboot:
    lines+=['The run paused for a reboot during CPU BiCGStab. The resumed '
            'process repeated setup and warmup, then reran the interrupted '
            'tile and continued from the accepted checkpoints. Solve totals '
            'sum the accepted tile times across both processes and include '
            'lid-close sleeps within measured tiles. The recorded sleeps add '
            'about 17 minutes to the 192-RHS CPU BiCGStab total and about 42 '
            'minutes to both CPU smoother totals. '
            '[Checkpoint and power log notes](RESUME_THREE_LEVEL.md).','']
lines+=[
        '| RHS | CPU BiCG / MG | CPU smoother / MG | Metal BiCG / MG | Metal smoother / MG | MG CPU / Metal |',
        '|---:|---:|---:|---:|---:|---:|']
for n in (96,192):
    ratios=[total(b,m,n)/total(b,'mg3',n) for b in processors for m in methods[:2]]+[total('cpu','mg3',n)/total('metal','mg3',n)]
    lines.append(f'| {n} | '+' | '.join(f'{v:.2f}×' for v in ratios)+' |')
lines+=['','## Convergence and measured work','',
       f'Fine FGMRES restart 8 and common inner target 0.002. MG uses 12 '
       'vectors per chirality at both interpolations, 4⁴ then 2⁴ aggregation, '
       'eight setup filtering sweeps at weight 0.65, and one pre/post smoothing '
       'step at weight 1 on each non-bottom level. The first coarse solve is '
       f'FGMRES(4), target 0.1, cap 8. The bottom is sparse BiCGStab, target '
       f'0.1, cap **{bottom}**. Smoother-only uses two fine clover-block Jacobi '
       'steps. BiCGStab retains the unscaled odd Schur baseline.','',
       '| RHS | Processor / method | Outer corrections/tile | Krylov steps/tile | Fine Schur calls (tile sum) | Level 1 SpMV (tile sum) | Bottom SpMV (tile sum) | Worst full FP64 residual |',
       '|---:|---|---:|---:|---:|---:|---:|---:|']
for n in (96,192):
    for b in processors:
        for m in methods:
            rs=[r for t in samples for r in selected(b,m,t,n)]
            def span(field):
                x=[int(r[field]) for r in rs]
                return str(min(x)) if min(x)==max(x) else f'{min(x)}–{max(x)}'
            work=[total(b,m,n,k) for k in ('fine_schur_calls','level1_products','level2_products')]
            lines.append(f'| {n} | {b} / {labels[m]} | {span("outer")} | {span("krylov")} | '+' | '.join(f'{v:g}' for v in work)+f' | {worst(b,m,n,"max_true_residual"):.3g} |')
lines+=['','Each Schur/SpMV call processes the resident tile. Sums at 192 RHS '
        f'contain {192//tile} tile workloads. Krylov steps sum per-RHS maxima across '
        'outer corrections; BiCGStab and FGMRES steps are different work. '
        'Metal bottom counts include padded/masked dispatches; CPU coarse '
        'iterations stop adaptively. Bottom residuals are diagnostic and '
        'never accept a full solve.','',
        'Full CPU FP64 residual history for the first tile (trial 0):','',
        '| Outer | CPU BiCG | CPU smoother | CPU MG | Metal BiCG | Metal smoother | Metal MG |',
        '|---:|---:|---:|---:|---:|---:|---:|']
history=defaultdict(dict)
for t in traces:
    if t['trial']=='0' and t['offset']=='0':history[int(t['outer'])][t['backend'],t['method']]=float(t['max_true_residual'])
for step,v in sorted(history.items()):lines.append(f'| {step} | '+' | '.join(f'{v[b,m]:.3g}' if (b,m) in v else 'accepted earlier' for b in processors for m in methods)+' |')
lines+=['','## Bottom-budget screening','',
        'The eight-RHS screening used numerical work and full residual '
        'acceptance to choose the bottom budget. Screening times varied '
        'with desktop load and are not final benchmark samples.','',
        '| Resident RHS | Bottom cap | Fine Krylov steps | Fine Schur | Bottom SpMV | Worst actual bottom residual | Full FP64 residual |',
        '|---:|---:|---:|---:|---:|---:|---:|']
for file in ('results/threelevel32-positive-pilot.csv','results/threelevel32-bottom8-pilot.csv'):
    with open(file) as stream:
        for r in csv.DictReader(stream):
            if r['backend']=='metal' and r['method']=='mg3':
                lines.append(f'| {r["tile_rhs"]} | {r["bottom_limit"]} | {r["krylov"]} | {r["fine_schur_calls"]} | {r["level2_products"]} | {float(r["worst_bottom_true_residual"]):.4g} | {float(r["max_true_residual"]):.3g} |')
lines+=['','The smaller cap preserves fine work while removing padded '
        'bottom products on that screening tile. The complete four-RHS '
        'sweeps above provide the final timing and convergence comparison.','',
        '## Time components and memory','',
        '| Processor / method | 192-RHS FP32 inner s | 192-RHS CPU FP64 residual s | Hierarchy setup s | Method workspace s | Process high-water RSS GiB | Largest recorded Metal allocation GiB |',
        '|---|---:|---:|---:|---:|---:|---:|']
for b in processors:
    for m in methods:
        rs=[r for t in samples for r in selected(b,m,t,192)]
        lines.append(f'| {b} / {labels[m]} | {total(b,m,192,"inner_ms")/1000:.2f} | {total(b,m,192,"residual_ms")/1000:.2f} | {st.median(float(r["hierarchy_ms"]) for r in rs)/1000:.2f} | {st.median(float(r["workspace_ms"]) for r in rs)/1000:.2f} | {max(float(r["peak_rss_mib"]) for r in rs)/1024:.2f} | {max(float(r["metal_allocated_mib"]) for r in rs)/1024:.2f} |')
lines+=['','RSS is the cumulative process high-water mark. Metal allocation '
        'counts the device resources currently allocated at each recorded '
        'tile; it is not a residency measurement. Shared memory appears in '
        'both accounting systems, so the columns must not be added. Workspace '
        'setup and hierarchy setup are reused across tiles. The common '
        'initial Metal library compilation is outside solve time.','',
        'The compact stencil avoids expanded CSR duplicates. The outer solver '
        'uses only b, x and residual as full FP64 fields, reusing the residual '
        'buffer for M64*x. Truth is regenerated for the independent solution '
        'check. Correction buffers are reused after their previous contents '
        'are no longer needed. FP64 hierarchy setup copies are discarded after '
        'validation; full FP64 gauge and clover data remain for all residual '
        'checks. See [THREE_LEVEL.md](THREE_LEVEL.md) for the precise layout.','',
        '## Mass location and harder pilot','',
        'The separate fixed-memory FP64 variational probe gives the following '
        'singular Rayleigh upper estimates. It recomputes the actual full '
        'Q² eigen-residual and passes free-field dispersion checks.','',
        '| Mass | Steps | Seed | Gap estimate | Actual Q² eigen-residual | Stationary mass | Reached 1e-5 eigen-residual? |',
        '|---:|---:|---:|---:|---:|---:|---|']
for f in ('results/compact-critical32.csv','results/compact-critical32-independent.csv'):
    with open(f) as stream:
        for s in csv.DictReader(stream):lines.append(f'| {float(s["mass"]):g} | {s["steps"]} | {s["seed"]} | {float(s["sigma_estimate"]):.7f} | {float(s["true_H2_residual"]):.3g} | {float(s["stationary_mass"]):.6f} | {"yes" if s["converged"]=="1" else "no"} |')
lines+=['','These approximate candidates locate a low-gap mass region; they '
        'do not reach the eigen-residual target and do not certify a global '
        'spectral minimum. The stationary mass minimizes ||M(m)v|| for that '
        'particular candidate. A physical PCAC critical mass remains '
        'uncalibrated. The full solve target remains strict regardless of '
        'these approximate spectral diagnostics.','',
        'At m=-0.12, an earlier four-RHS Metal MG pilot passed in 30 Krylov '
        'steps and 94 fine Schur calls, with full residual 8.4e-13. Plain '
        'Metal BiCGStab did not meet its inner target within 2,048 iterations '
        'in that pilot; it never returned an accepted full solution. The '
        f'complete comparison above uses m={mass:g} on the more positive '
        'side of the screened region. Pilot times are not included in the '
        'complete comparison.','',
        '## Validation and reproduction','',
        f'All **{len(rows)} measured tiles**, or **{len(rows)*tile:,} RHS '
        f'executions**, passed full original-equation FP64 acceptance. '
        f'Worst final residual {max(float(r["max_true_residual"]) for r in rows):.4g}; '
        f'worst independent known-solution error {max(float(r["solution_error"]) for r in rows):.4g}. '
        'Every final CSV residual was checked against its per-outer trace. '
        'The 32⁴ projections pass P†P=I/isometry, direct RSP/RAP and gamma5 '
        'tests at about 2e-15, and CPU/Metal restriction, prolongation and '
        'sparse products agree within 1e-6. Small-lattice tests compare the '
        'compact full and Schur operators against the original implementation, '
        'exercise all six solver paths and zero RHS, and check resumption. '
        'The original two-level solver passes its regression tests.','',
        '```sh','make build/metal_three_level build/compact_probe',
        'VECLIB_MAXIMUM_THREADS=1 caffeinate -i ./build/metal_three_level \\',
        f'  --size 32 --mass {mass:g} --csw 1 --tile-rhs {tile} --rhs 96,192 --samples {len(samples)} \\',
        f'  --threads 8 --coarse-iters 8 --bottom-iters {bottom} \\',
        '  --csv results/threelevel32-compare.csv',
        'python3 threelevel_summary.py'+(' --resumed-after-reboot' if args.resumed_after_reboot else ''),'```','',
        'Add `--resume` to the identical benchmark command to skip accepted '
        'tiles. All code is confined to the standalone prototype. Chroma and '
        'SuperbBLAS sources are unchanged.','',
        f'Raw [measurements]({args.csv.as_posix()}), '
        f'[full FP64 traces]({args.csv.as_posix()}.trace.csv), '
        '[run log](results/threelevel32-compare.log), '
        '[methodology and source layout](THREE_LEVEL.md).','']
args.output.write_text('\n'.join(lines))
print(f'Wrote {args.output}: audited {len(rows)} tiles, {len(samples)} sweep(s)')
