#!/usr/bin/env python3
"""Audit and summarize the negative-mass CPU/Metal comparison."""
import argparse
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

ap=argparse.ArgumentParser()
ap.add_argument('csv',nargs='?',type=Path,default=Path('results/nearcritical-compare16.csv'))
ap.add_argument('--output',type=Path,default=Path('NEAR_CRITICAL_RESULTS.md'))
args=ap.parse_args()
with args.csv.open() as f:rows=list(csv.DictReader(f))
with Path(str(args.csv)+'.trace.csv').open() as f:traces=list(csv.DictReader(f))
groups=defaultdict(list)
for r in rows:
    assert int(r['extent'])==16 and float(r['csw'])==1 and float(r['mass'])<0
    assert int(r['tile_rhs'])==96 and int(r['rhs'])==96*int(r['tiles'])
    assert int(r['coarse_limit'])==16 and float(r['coarse_tolerance'])==.1
    assert int(r['bi_chunk'])==1 and int(r['restart'])==8 and int(r['vectors'])==12
    assert int(r['workers'])==8 and int(r['coarse_dim'])==6144 and int(r['coarse_nnz_blocks'])==8448
    assert float(r['inner_tolerance'])==.002 and float(r['tolerance'])==1e-10
    assert all(math.isfinite(float(r[k])) for k in ('total_ms','max_true_residual','solution_error','worst_coarse_true_residual'))
    assert float(r['max_true_residual'])<=1e-10 and float(r['solution_error'])<=1e-7
    if not r['backend'].endswith('sparse_mg'):assert int(r['coarse_spmv_calls'])==0 and float(r['hierarchy_setup_ms'])==0
    groups[int(r['rhs']),r['backend']].append(r)
masses={float(r['mass']) for r in rows};assert len(masses)==1
mass=next(iter(masses));assert mass==-.12
methods=('bicgstab','smooth','sparse_mg')
backends=tuple(p+'_'+m for p in ('cpu','metal') for m in methods)
sample_counts={len(groups[rhs,b]) for rhs in (96,192) for b in backends}
assert len(sample_counts)==1
samples=next(iter(sample_counts));assert samples in (1,3,5,7,9)
assert len(rows)==12*samples
last=defaultdict(dict)
for t in traces:
    assert int(t['extent'])==16 and float(t['mass'])==mass
    assert math.isfinite(float(t['max_true_residual']))
    key=t['backend'],int(t['trial']);off=int(t['offset'])
    if off not in last[key] or int(t['outer_step'])>int(last[key][off]['outer_step']):last[key][off]=t
for r in rows:
    selected=[t for off,t in last[r['backend'],int(r['trial'])].items() if off<int(r['rhs'])]
    assert len(selected)==int(r['tiles'])
    assert all(float(t['max_true_residual'])<=1e-10 for t in selected)
    assert math.isclose(max(float(t['max_true_residual']) for t in selected),float(r['max_true_residual']),rel_tol=1e-10)
physical=sum(len(v) for v in last.values());assert physical==12*samples

def vals(rhs,b,field):return [float(r[field]) for r in groups[rhs,b]]
def med(rhs,b,field='total_ms'):return st.median(vals(rhs,b,field))
def span(rhs,b,field):
    v=sorted({int(r[field]) for r in groups[rhs,b]})
    return str(v[0]) if len(v)==1 else f'{v[0]}–{v[-1]}'
labels={'bicgstab':'BiCGStab','smooth':'smoother FGMRES','sparse_mg':'sparse MG'}
with Path('results/critical-refine.csv').open() as f:spectrum=list(csv.DictReader(f))
with Path('results/critical-independent.csv').open() as f:independent=list(csv.DictReader(f))
with Path('results/solver-compare16.csv').open() as f:reference=list(csv.DictReader(f))
ref=defaultdict(list)
for r in reference:
    if float(r['mass'])==.1:ref[int(r['rhs']),r['backend']].append(r)
lines=['# Clover-Wilson comparison near the estimated gap minimum','',
       f'Tested bare mass **m₀={mass:g}** on the same 16⁴ synthetic gauge field '
       '(strength 0.25), c_sw=1, QDP conventions, and antiperiodic time. '
       'The chosen mass is on the positive side of a finite-configuration '
       'low-gap minimum estimated near -0.14. Full residuals of the original '
       'reconstructed clover-Wilson equation are computed on the CPU in FP64 '
       'for every component and every RHS, with target 1e-10. Inner work '
       'and even/odd preparation/reconstruction are FP32 on the selected '
       'processor. Hardware: Apple M1 Max, 24 GPU cores, 32 GB, macOS 26.7.1. '
       'Built with Apple clang 21.0.0 using `-O3 -mcpu=apple-m1 '
       '-ffp-contract=fast`; full flags are in the Makefile.','',
       '## Complete solve times','',
       (f'Medians of {samples} paired trials after a full resident-tile warmup per method. '
        if samples>1 else 'One paired timing trial after a full resident-tile warmup per method. ')+
       '**96 RHS are resident; 192 RHS use two sequential 96-RHS tiles.** '
       'The first actual measured tile supplies the 96-RHS result and the '
       'sum of both distinct measured tiles supplies the 192-RHS result. '
       'Every tile starts from zero. Metal runs first, followed by CPU; '
       'method workspaces are allocated sequentially, retaining the same '
       'operator, hierarchy and common FP64 outer buffers.','',
       '| RHS | CPU BiCG s | CPU smoother FG s | CPU MG s | Metal BiCG s | Metal smoother FG s | Metal MG s |',
       '|---:|---:|---:|---:|---:|---:|---:|']
for rhs in (96,192):lines.append(f'| {rhs} | '+' | '.join(f'{med(rhs,b)/1000:.3f}' for b in backends)+' |')
lines+=['','Solve time includes resets, all fine/coarse products, vector '
        'operations/reductions, smoothing, small CPU QR, GPU submission/waits, '
        'FP32 preparation and reconstruction, FP64 updates, and full CPU '
        'FP64 residual checks. Hierarchy/workspace/problem setup and '
        'manufactured RHS construction are separate. Initial common Metal '
        'library compilation is outside solve timing. Known-solution error '
        'checking follows residual acceptance outside timing. CPU uses '
        'eight workers; setup BLAS is single-threaded. These timings have '
        'ordinary local system activity and no thermal/scheduling isolation.','',
        '| RHS | CPU BiCG / MG | CPU smoother / MG | Metal BiCG / MG | Metal smoother / MG | MG CPU / Metal |',
        '|---:|---:|---:|---:|---:|---:|']
for rhs in (96,192):
    gains=[med(rhs,p+'_'+m)/med(rhs,p+'_sparse_mg') for p in ('cpu','metal') for m in ('bicgstab','smooth')]
    gains.append(med(rhs,'cpu_sparse_mg')/med(rhs,'metal_sparse_mg'))
    lines.append(f'| {rhs} | '+' | '.join(f'{v:.2f}×' for v in gains)+' |')
lines+=['','The first four ratios compare each baseline with MG; the last '
        'is the Metal speedup over CPU for MG. Medians and timing ranges describe '
        'these runs; small differences within the ranges need caution.','',
        '| RHS | Method | Min–max solve s | Full CPU FP64 residual s | FP32 inner s |',
        '|---:|---|---:|---:|---:|']
for rhs in (96,192):
    for b in backends:
        v=vals(rhs,b,'total_ms')
        lines.append(f'| {rhs} | {b} | {min(v)/1000:.3f}–{max(v)/1000:.3f} | {med(rhs,b,"residual_ms")/1000:.3f} | {med(rhs,b,"inner_ms")/1000:.3f} |')
lines+=['',f'For Metal MG at 96 RHS, full CPU FP64 residual checks take '
        f'{med(96,"metal_sparse_mg","residual_ms")/1000:.3f} s, or '
        f'{100*med(96,"metal_sparse_mg","residual_ms")/med(96,"metal_sparse_mg"):.1f}% '
        'of complete solve time. This cost remains in every timing above.','',
        '## Convergence and actual work','',
        'The inner target is 0.002 for every method. Fine FGMRES restart '
        'remains 8, and smoother-only uses two clover-block Jacobi steps at '
        'weight 1.0. MG has one pre/post step at the same weight, 12 vectors '
        'per chirality, 4⁴ aggregates, and eight setup filtering sweeps '
        'at weight 0.65. Its sparse coarse solve uses a stronger 16-iteration '
        'cap and target preconditioned residual 0.1, chosen by the pilot '
        'below. BiCGStab retains its odd-Schur baseline, and Metal checks '
        'convergence every iteration (chunk=1). All methods have the '
        'same full FP64 target.','',
        '| RHS | Method | Outer corrections/tile | Summed Krylov steps/tile | Total fine S calls | Coarse SpMV | Worst full FP64 residual |',
        '|---:|---|---:|---:|---:|---:|---:|']
for rhs in (96,192):
    for b in backends:
        omin,omax=min(vals(rhs,b,'outer_min')),max(vals(rhs,b,'outer_max'))
        kmin,kmax=min(vals(rhs,b,'krylov_min')),max(vals(rhs,b,'krylov_max'))
        outer=f'{omin:g}' if omin==omax else f'{omin:g}–{omax:g}'
        kry=f'{kmin:g}' if kmin==kmax else f'{kmin:g}–{kmax:g}'
        lines.append(f'| {rhs} | {b} | {outer} | {kry} | {span(rhs,b,"fine_schur_calls")} | {span(rhs,b,"coarse_spmv_calls")} | {max(vals(rhs,b,"max_true_residual")):.3g} |')
lines+=['','Krylov counts sum per-RHS maxima across outer corrections. '
        'BiCGStab and FGMRES iterations are different work. Every fine '
        'Schur application is counted, including smoothing, refreshes '
        'and GPU half-step dispatches. Calls at 192 RHS sum both 96-wide '
        'tiles. GPU coarse counts include masked/padded dispatches. '
        'Actual coarse residuals are diagnostic and do not accept a solve.','',
        'Example full FP64 residual history for global RHS offset 96, trial 0:','',
        '| Outer correction | CPU BiCG | CPU smoother FG | CPU MG | Metal BiCG | Metal smoother FG | Metal MG |',
        '|---:|---:|---:|---:|---:|---:|---:|']
history=defaultdict(dict)
for t in traces:
    if t['trial']=='0' and t['offset']=='96':history[int(t['outer_step'])][t['backend']]=float(t['max_true_residual'])
for step,v in sorted(history.items()):lines.append(f'| {step} | '+' | '.join(f'{v[b]:.3g}' if b in v else 'accepted earlier' for b in backends)+' |')
lines+=['','### Comparison with the earlier m₀=+0.1 test','',
        'Historical positive-mass data provide a convergence/work reference. '
        'These are not additional fresh timing samples. The near-critical '
        'MG configuration also has a stronger coarse budget; that change '
        'is part of the reported configuration.','',
        '| Method | Old +0.1 fine S / 96 RHS | New -0.12 fine S / 96 RHS | Old +0.1 Krylov steps | New -0.12 Krylov steps |',
        '|---|---:|---:|---:|---:|']
for b in backends:
    old_s=sorted({int(r['fine_schur_calls']) for r in ref[96,b]});old_k=sorted({int(r['krylov_sum']) for r in ref[96,b]})
    lines.append(f'| {b} | {"/".join(map(str,old_s))} | {span(96,b,"fine_schur_calls")} | {"/".join(map(str,old_k))} | {span(96,b,"krylov_sum")} |')
lines+=['','## Spectral location and its limits','',
        'The separate FP64 diagnostic applies `Q=gamma5*M` and estimates '
        'the low spectrum of `Q²=M†M` using fully reorthogonalized Lanczos. '
        'The gap convention is the minimum absolute Q eigenvalue '
        '([Del Debbio et al., section 2](https://arxiv.org/abs/hep-lat/0512021)). '
        'Candidate vectors are reconstructed and their actual full-lattice '
        'FP64 Q² eigen-residual is evaluated. Free-field tests match the '
        'analytic Wilson dispersion.','',
        '| Bare mass | Steps | Start seed | Singular Rayleigh estimate | Actual Q² eigen-residual | Stationary-mass diagnostic | Reached 1e-5 eigen-residual? |',
        '|---:|---:|---:|---:|---:|---:|---|']
for s in spectrum+independent:
    lines.append(f'| {float(s["mass"]):g} | {s["steps"]} | {s["seed"]} | {float(s["sigma_min_estimate"]):.7f} | {float(s["true_H2_residual"]):.3g} | {float(s["stationary_mass"]):.6f} | {"yes" if s["converged"]=="1" else "no"} |')
lines+=['','The scan estimates a low-gap minimum near -0.14, with m=-0.12 '
        'on its positive side. Different starts/lengths near -0.14 give '
        'consistent approximate Rayleigh values, but the finite-volume '
        'runs do not reach their requested eigen-residual tolerance. '
        'Their values are Rayleigh upper estimates, so no precise '
        'critical mass or certified global spectral bound is claimed. '
        'The stationary-mass diagnostic is `m-Re<v,M(m)v>` for the '
        'normalized candidate; it locates the minimum of the mass '
        'quadratic for that candidate. This is a fixed-field spectral '
        'proxy, with antiperiodic time retained. A physical PCAC critical '
        'mass has not been determined. The finite synthetic field still '
        'has a substantial estimated gap near the minimum. **The full '
        'solve acceptance remains 1e-10 in FP64 regardless of these '
        'approximate spectral diagnostics.**','',
        '## Coarse solver screening','',
        'The four-step coarse cap from the positive-mass test was inadequate '
        'near this mass. The small-RHS screening uses the same field and '
        'strict full residual checks; it selects coarse accuracy and '
        'fine work, not a final timing result. Some screening diagnostics '
        'overlapped other functional/spectral checks. The final comparison '
        'runs methods sequentially.','',
        '| Coarse cap | Resident RHS | Fine Krylov steps | Fine S | Worst actual coarse residual | Full FP64 residual |',
        '|---:|---:|---:|---:|---:|---:|']
with Path('results/nearcritical-pilot.csv').open() as f:pilots=[r for r in csv.DictReader(f) if float(r['mass'])==mass and r['backend']=='metal_sparse_mg']
for cap in (8,12,16):
    with Path(f'results/nearcritical-coarse{cap}.csv').open() as f:pilots+=list(csv.DictReader(f))
for r in pilots:
    cap=r.get('coarse_limit','4')
    lines.append(f'| {cap} | {r["tile_rhs"]} | {r["krylov_sum"]} | {r["fine_schur_calls"]} | {float(r["worst_coarse_true_residual"]):.4g} | {float(r["max_true_residual"]):.3g} |')
lines+=['','Caps 12 and 16 gave the same pilot convergence with a good '
        'actual coarse residual; 16 was chosen for the wider RHS batch. '
        'CPU stops submitting coarse products after convergence. Metal '
        'encodes the cap with per-RHS masking, so extra padded work is '
        'inside the measured MG time. FGMRES stores the preconditioned '
        'basis to permit the inexact/variable correction.','',
        '## Setup and many RHS','',
        '| RHS | MG hierarchy s | CPU hierarchy-only break-even vs best baseline | Metal hierarchy-only break-even vs best baseline |',
        '|---:|---:|---:|---:|']
for rhs in (96,192):
    setup=med(rhs,'cpu_sparse_mg','hierarchy_setup_ms');be=[]
    for p in ('cpu','metal'):
        best=min(med(rhs,p+'_'+m) for m in ('bicgstab','smooth'));saved=best-med(rhs,p+'_sparse_mg')
        be.append(f'{setup/saved:.2f} batches' if saved>0 else 'No warm gain')
    lines.append(f'| {rhs} | {setup/1000:.3f} | {be[0]} | {be[1]} |')
lines+=['','These estimates divide hierarchy setup by observed warm batch '
        'time saved, excluding differences in method-specific workspace '
        'setup and first-use costs. The hierarchy is reused for all '
        'RHS, processors and trials at this mass; more RHS amortize it.','',
        '| Method | Inner workspace s | Common outer buffers s | Common problem setup s | Peak process RSS GiB | Worst coarse true residual |',
        '|---|---:|---:|---:|---:|---:|']
for b in backends:
    lines.append(f'| {b} | {med(96,b,"inner_workspace_setup_ms")/1000:.3f} | {med(96,b,"outer_workspace_setup_ms")/1000:.3f} | {med(96,b,"problem_setup_ms")/1000:.3f} | {max(vals(192,b,"peak_rss_mib"))/1024:.2f} | {max(vals(192,b,"worst_coarse_true_residual")):.3g} |')
lines+=['','RSS is the cumulative process high-water mark, not an isolated '
        'per-method measurement or total system memory. The sparse '
        'matrix retains 6,144 coarse unknowns and 37.125 MiB of FP32 '
        'coefficients; large fine vectors dominate memory.','',
        '## Validation and reproduction','',
        f'All {len(rows)} paired logical records passed every-RHS full CPU '
        f'FP64 residual acceptance. They derive from {physical} physical '
        f'96-wide tile solves ({physical*96:,} RHS executions). Worst final '
        f'relative residual: {max(float(r["max_true_residual"]) for r in rows):.3g}; '
        f'maximum known-solution error: {max(float(r["solution_error"]) for r in rows):.3g}. '
        'CSV acceptance values were audited against final per-tile traces. '
        'Negative-mass small-lattice tests cover all six paths and zero '
        'RHS. The actual 16⁴ m=-0.12 sparse hierarchy passes full FP64 '
        'Galerkin/gamma5 checks (relative Galerkin error about 1.8e-15), '
        'plus strict CPU/Metal coarse solves. Local clover Cholesky '
        'factors are positive; the reported smallest even-factor pivot '
        'is about 3.4854.','',
        '```sh','make build/critical_probe build/metal_solver_compare',
        'VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \\',
        f'  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass -.12 --samples {samples} \\',
        '  --threads 8 --csw 1 --bi-chunk 1 --coarse-iters 16 --coarse-tol .1 \\',
        '  --metal-first --csv results/nearcritical-compare16.csv',
        'python3 nearcritical_summary.py','```','',
        'See [NEAR_CRITICAL.md](NEAR_CRITICAL.md) for spectral methodology, '
        'pilot details and the timing protocol. These are synthetic-field '
        'results. Chroma and SuperbBLAS sources are unchanged.','',
        f'Raw [measurements]({args.csv.as_posix()}), '
        f'[full FP64 traces]({args.csv.as_posix()}.trace.csv), '
        '[run log](results/nearcritical-compare16.log), '
        '[spectral scan](results/critical-refine.csv), '
        '[independent spectral check](results/critical-independent.csv).','']
args.output.write_text('\n'.join(lines))
print(f'Wrote {args.output}: {len(rows)} records, {physical} tile solves, {samples} trials')
