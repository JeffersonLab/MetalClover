#!/usr/bin/env python3
"""Render the accepted real 32-cubed result; all rates have explicit scope."""
import json
import math
from pathlib import Path
import sys

root = Path(sys.argv[1])
s = json.loads((root / 'summary.json').read_text())
mg = s['methods']['near-null-svd']
bi = s['methods']['metal-bicgstab']
old = 916.269253
old_setup = 269.226134
rows = []
for name, item in [('Near-null/SVD Metal MG', mg), ('Metal BiCGStab', bi)]:
    rows.append(f"| {name} | {item['mean_seconds_per_rhs']:.3f} | "
                f"{sum(item['fine_schur_applications']) / 2:.1f} | "
                f"{item['weighted_operator_GFLOP_per_second']:.2f} | "
                f"{item['maximum_original_Chroma_residual']:.4e} |")
energy = s['basis_svd']
text = f'''# Near-null/SVD MG on the real 32³×256 light-mass field

The **mass −0.0860** configuration solves on the **32 GiB Apple M1 Max**
with two deterministic complex Z4 RHS (tags 20261107 and 20261108), one
resident GPU RHS at a time. One MG hierarchy was built and reused for both.
Every field passed the independent original Chroma CPU FP64 full-system
check at **1e-12**. This measures **two RHS**, not 96 or 192.

| Method | Mean solve (s/RHS) | Fine Schur applications/RHS | Nominal operator GFLOP/s | Maximum original Chroma FP64 residual |
|---|---:|---:|---:|---:|
{chr(10).join(rows)}

The solve time includes GPU reduction/reconstruction, FP64 solution updates,
slab writeback and full CPU FP64 residual checks. Import, common operator
setup/probe and subsequent independent native verification are outside this
solve time. The rate is a **nominal source-level dispatched-operator model**,
not hardware-counter throughput or total algorithm FLOP/s. Its numerator
excludes Krylov-vector/reduction and CPU FP64 operations; their time remains
in the denominator. Padded coarse work is counted. See
[OPERATOR_WORK.md](OPERATOR_WORK.md) for formulas and interfaces.

MG is **{s['warm_speedup']:.2f}× faster** on solve time, using
**{100*s['fine_work_reduction']:.2f}% fewer fine Schur applications**.
MG took {mg['seconds'][0]:.3f}/{mg['seconds'][1]:.3f} seconds and
{mg['krylov_iterations'][0]}/{mg['krylov_iterations'][1]} fine FGMRES iterations;
BiCGStab took {bi['seconds'][0]:.3f}/{bi['seconds'][1]:.3f} seconds and
{bi['krylov_iterations'][0]}/{bi['krylov_iterations'][1]} iterations.
Each solve used five defect corrections. The independent native residuals
rather than a coarse, recursive or FP32 residual determine acceptance.

## Setup and many RHS

MG setup cost **{mg['hierarchy_seconds']:.3f} seconds**, including 4,272 fine
Schur applications. Including this setup amortized over the two measured
RHS gives **{mg['measured_seconds_per_rhs_including_setup']:.3f} s/RHS**,
which is **{s['speedup_including_measured_setup']:.2f}× faster** than BiCGStab.
The following are arithmetic projections of that same measured hierarchy
and mean solve time; additional RHS were not timed:

| Reused RHS | Projected MG seconds/RHS including hierarchy |
|---:|---:|
| 96 | {mg['projected_seconds_per_rhs_including_setup']['96']:.3f} |
| 192 | {mg['projected_seconds_per_rhs_including_setup']['192']:.3f} |

The earlier filtered-basis result on this same field was {old:.3f} s/RHS
with {old_setup:.3f} s setup ([historical result](LOW_MEMORY.md)). The new
warm solve is {old/mg['mean_seconds_per_rhs']:.2f}× faster than that historical
MG run. Its additional setup amortizes after approximately
{(mg['hierarchy_seconds']-old_setup)/(old-mg['mean_seconds_per_rhs']):.2f} RHS.
The filtered-basis solver was not retimed in this experiment.

## Bounded basis setup

Fine aggregates are 4×4×4×4, with eight vectors per chirality retained at
both transfers. The fine restart is eight. The fine candidates use 32 local
warm steps, homogeneous FGMRES restart eight with a 128-iteration cap and
0.001 target. All 24 candidates hit that cap; their achieved reductions are
recorded individually. The second-level filter uses 64 steps; the bottom
sparse iteration cap is eight. No dense coarse solve is used.

At this volume a single FP32 odd candidate occupies 384 MiB. Holding all
24 candidates would require 9 GiB. The incremental construction instead
uses a **13-column pool (4.875 GiB)**: compress to nine weighted left singular
vectors `U*Sigma`, append four new candidates, and repeat. The final local
FP64 SVD retains eight unweighted left singular vectors. Carrying singular
values preserves the represented covariance weights between batches.
The fine setup pool is reduced by **4.125 GiB**.

The final CPU FP64 local basis has maximum `|U†U−I|` of
**{energy['orthogonality']:.4e}**, before casting for the GPU. Its energy
lower bound relative to the original 24 normalized candidates is
**{100*energy['energy_lower_bound']:.3f}%**, up to FP32 carry/normalization
roundoff. This is a lower bound from the represented covariance, not the
retained-energy fraction of a full 24-column SVD: previously discarded
positive covariance contributions can be recaptured by later directions.
The final truncation retains {100*energy['energy']:.3f}% of the represented
energy; the minimum local fraction is {100*energy['minimum_energy']:.3f}%.

Measured Metal allocation at ready is **{mg['GPU_allocation_at_ready_GiB']:.3f}
GiB** for MG and **{bi['GPU_allocation_at_ready_GiB']:.3f} GiB** for BiCGStab,
against a 24.960 GiB device recommendation. MG process peak RSS is
{mg['peak_RSS_GiB']:.3f} GiB. These unified-memory metrics overlap and must
not be added. CPU FP64 clover coefficients are regenerated from original
FP64 gauge links; the solution is file backed and residual work is streamed
in eight-time-slice slabs. GPU inversion and original native verification
ran sequentially. This result establishes the standalone solver footprint;
Chroma also owns its gauge/clover/spinor storage.

## Operator and evidence

The header action uses xi0=4.3, nu=1.265, clovCoeffR=1.58932722549812 and
clovCoeffT=0.902783591772098, two spatial stout steps at rho=0.14, and
antiperiodic time. No extra clover rescaling is introduced. The whole-lattice
CPU FP64 operator probe agrees with exported original Chroma to 1.7086e-16;
the native checker reproduces that export exactly.

[Structured results](results/near-null-svd32/summary.json),
[MG timings](results/near-null-svd32/near.mg.csv),
[BiCGStab timings](results/near-null-svd32/bicg.bicg.csv),
[MG native checks](results/near-null-svd32/mg.native.csv), and
[BiCGStab native checks](results/near-null-svd32/bicg.native.csv) contain the
accepted measurements. Logs, full refinement traces, resource samples,
action XML, source/binary hashes and the bounded-SVD diagnostics are stored
beside them. Large binary fields remain generated local scratch artifacts.

The immutable numerical benchmark binary predates the installable API and
namespace refactor. Its exact source is recoverable from base commit
369c3dc plus `benchmark-source.patch.gz` (decompressed); the source and binary hashes are
recorded. The final namespaced driver is checked on the real 12³×96 field
against original Chroma, including comparison of saved solution bits with
the pre-refactor control. The installable C ABI and native Chroma array
harness have separate real-field checks.
'''
Path(sys.argv[2]).write_text(text)
