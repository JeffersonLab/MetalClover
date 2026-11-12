# Near-null/SVD MG on the real 32³×256 light-mass field

The **mass −0.0860** configuration solves on the **32 GiB Apple M1 Max**
with two deterministic complex Z4 RHS (tags 20261107 and 20261108), one
resident GPU RHS at a time. One MG hierarchy was built and reused for both.
Every field passed the independent original Chroma CPU FP64 full-system
check at **1e-12**. This measures **two RHS**, not 96 or 192.

| Method | Mean solve (s/RHS) | Fine Schur applications/RHS | Nominal operator GFLOP/s | Maximum original Chroma FP64 residual |
|---|---:|---:|---:|---:|
| Near-null/SVD Metal MG | 474.499 | 648.0 | 131.85 | 1.7355e-14 |
| Metal BiCGStab | 2189.472 | 12714.0 | 381.28 | 2.5931e-14 |

The solve time includes GPU reduction/reconstruction, FP64 solution updates,
slab writeback and full CPU FP64 residual checks. Import, common operator
setup/probe and subsequent independent native verification are outside this
solve time. The rate is a **nominal source-level dispatched-operator model**,
not hardware-counter throughput or total algorithm FLOP/s. Its numerator
excludes Krylov-vector/reduction and CPU FP64 operations; their time remains
in the denominator. Padded coarse work is counted. See
[OPERATOR_WORK.md](OPERATOR_WORK.md) for formulas and interfaces.

MG is **4.61× faster** on solve time, using
**94.90% fewer fine Schur applications**.
MG took 470.209/478.788 seconds and
205/209 fine FGMRES iterations;
BiCGStab took 2154.528/2224.416 seconds and
6263/6451 iterations.
Each solve used five defect corrections. The independent native residuals
rather than a coarse, recursive or FP32 residual determine acceptance.

## Setup and many RHS

MG setup cost **964.986 seconds**, including 4,272 fine
Schur applications. Including this setup amortized over the two measured
RHS gives **956.992 s/RHS**,
which is **2.29× faster** than BiCGStab.
The following are arithmetic projections of that same measured hierarchy
and mean solve time; additional RHS were not timed:

| Reused RHS | Projected MG seconds/RHS including hierarchy |
|---:|---:|
| 96 | 484.551 |
| 192 | 479.525 |

The earlier filtered-basis result on this same field was 916.269 s/RHS
with 269.226 s setup ([historical result](LOW_MEMORY.md)). The new
warm solve is 1.93× faster than that historical
MG run. Its additional setup amortizes after approximately
1.57 RHS.
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
**8.2157e-15**, before casting for the GPU. Its energy
lower bound relative to the original 24 normalized candidates is
**96.555%**, up to FP32 carry/normalization
roundoff. This is a lower bound from the represented covariance, not the
retained-energy fraction of a full 24-column SVD: previously discarded
positive covariance contributions can be recaptured by later directions.
The final truncation retains 98.614% of the represented
energy; the minimum local fraction is 95.355%.

Measured Metal allocation at ready is **23.075
GiB** for MG and **9.285 GiB** for BiCGStab,
against a 24.960 GiB device recommendation. MG process peak RSS is
15.929 GiB. These unified-memory metrics overlap and must
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

The native Chroma adapter at commit `0c73a0c05d6e127898dbea1841ca296cbd45c8e1` on
`rge/devel-datalib-reorg` built from that branch's requested base. The real
12³×96 harness passed all 20 scalar, reconstruction, qprop-array and actual
SuperbBLAS colorvector/spin checks; its maximum full native FP64 residual was
2.3534e-14. See [native integration evidence](results/near-null-svd32/chroma-integration/validation.json).
This is separate from the standalone 32³×256 timing above.
