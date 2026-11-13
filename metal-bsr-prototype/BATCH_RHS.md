# Batched RHS support for real Chroma clover inversion

The Chroma `qprop` array interface now submits **genuine GPU RHS tiles** to
`MetalBiCGClover` and `MetalMGClover`. Fine Schur products, smoothers,
restriction/prolongation, sparse coarse products and bottom iterations carry
the RHS dimension. The inverters use independent Krylov solves/convergence
per column, sharing one gauge field, packed LDL clover data and MG hierarchy.
This does not introduce a coupled block Krylov algorithm.

These timings used the original independent full-field random spinors.
They count physical spinor RHS. The current distillation harness instead
uses four source spins per time-localized color vector and takes
`BASE_NUM_VECS`; 24 vectors mean 96 RHS and 48 mean 192. See
[DISTILLATION_RHS.md](DISTILLATION_RHS.md). The results below retain their
original source definition and have not been relabeled as distillation timings.

The real **12³×96**, mass **−0.0743** LIME field was tested on the **32 GiB
Apple M1 Max with 24 GPU cores**, one MPI rank, MetalClover CPU workers=8,
`OMP_NUM_THREADS=1`, `VECLIB_MAXIMUM_THREADS=1`. The header action uses
xi0=4.3, nu=1.265, clovCoeffR=1.58954576550983,
clovCoeffT=0.902866347885228, two spatial stout steps at rho=0.14 and
antiperiodic time. Every physical source is repeatable complex Z4 with seed
tag `20261113 + column` and unit full-field norm.

## Measured throughput

All three modes below solve the **same 96 sources**, with zero initial
solutions. Each method retains the **same qprop object**, and MG retains its
**one hierarchy**, across all modes and the subsequent 192-source run.
`GPU RHS` is the actual observed width, not merely an array API length.

| Method | GPU RHS | Complete qprop s/RHS | Speedup over scalar | Nominal operator GFLOP/s |
|---|---:|---:|---:|---:|
| Metal BiCGStab | 1 | 3.358 | 1.00× | 281.6 |
| Metal BiCGStab | 8 | 2.295 | 1.46× | 465.6 |
| Metal BiCGStab | 16 | 2.201 | 1.53× | 502.5 |
| Metal MG | 1 | 2.488 | 1.00× | 161.8 |
| Metal MG | 8 | 1.416 | 1.76× | 293.7 |
| Metal MG | 16 | 1.382 | 1.80× | 301.7 |

The additional **192-RHS** measurements use tiles of sixteen; a 192-RHS
scalar baseline was not timed:

| Method | GPU RHS | Complete qprop s/RHS | Nominal operator GFLOP/s |
|---|---:|---:|---:|
| Metal BiCGStab | 16 | 2.214 | 502.5 |
| Metal MG | 16 | 1.395 | 298.9 |

Native qprop construction took **19.714 s** for MG and
**1.175 s** for BiCGStab. It is outside solve timing.
One MG hierarchy serves all four measurements: 480 inversions over 192
unique sources. Common LIME input/stout-state creation, source generation,
reference-file I/O and additional independent post-timing verification are
outside solve timing. Native qprop RHS reduction, CPU FP64 refinement checks,
reconstruction and full residual evaluation are inside. First-use batch
workspace/pipeline creation is also inside and separately logged.

These are single wall-clock measurements, rather than repeated-run medians.
For the same 96 sources, all three widths had identical per-source iteration
counts: BiCGStab required 336–400 iterations and MG 76–78. At 192 sources,
the ranges were 336–408 and 76–78 respectively. Iterations are accumulated
over FP64 refinement corrections, with the same total budget per column.

Every source passed **full original Chroma CPU FP64 relative residual <=1e-12**.
The maximum across the measurements was **2.4328e-14**. Complete solutions
for the first 96 sources were compared against scalar BiCGStab, including
MG/scalar/batch agreement; maximum relative solution difference was
**6.3424e-13**. All 192 sources received the full native check.
Peak process RSS was **14.547 GiB**. The run completed
under a 26 GiB RSS limit and 4 GiB disk-reserve guard. The large scalar
reference file is a generated scratch artifact, excluded from Git.

The rate is **nominal dispatched GPU operator work / complete qprop wall
time**. FMA counts as two and a complex multiply/add as eight. It includes
inactive/padded lanes and the sparse coarse work. CPU FP64 and Krylov-vector
arithmetic are excluded from the numerator; their elapsed time is included.
These are not hardware-counter or total-algorithm FLOP rates. See
[OPERATOR_WORK.md](OPERATOR_WORK.md) for formulas and attribution.

## Interfaces and controls

Chroma keeps its existing `SystemSolver<LatticeFermion>` vector overload,
with solution/source vectors of shared pointers and one result per RHS.
The asymmetric even/odd qprop wrapper reduces sources, forwards the whole
array, reconstructs full solutions and evaluates native full residuals.
The Metal adapter calls the array C ABI once rather than looping over scalar
GPU solves. The existing SuperbBLAS `SB::ChimeraSolver`/`SB::doInversion`
path reaches the same implementation.

`<RHSBatch>16</RHSBatch>` sets the device cap; the default is 8, and 1 selects
scalar execution. A conservative device-allocation budget can lower the
actual width, retaining the large-volume scalar path. Gauge/clover, transfer
weights and sparse matrices share allocations across RHS widths; changing
width allocates workspaces without rebuilding the hierarchy. The caller's
`max_rhs` controls the group reaching qprop and its CPU temporary array.
Use bounded caller groups rather than retaining every full spinor at once.
MPI decomposition and concurrent calls on one context remain unsupported.

The C ABI remains version 1, with unchanged existing options/info layouts.
Its additive batch extension provides `metalclover_set_rhs_batch` and
`metalclover_get_batch_info`, advertised by `METALCLOVER_BATCH_API_VERSION`.
Returned column times are amortized batch shares; aggregate work/rate is
preferred for throughput. Per-column guesses, norms, residuals, iteration
counts and total budgets are independent. Residuals are normalized before
FP32 conversion and corrections rescaled in FP64, supporting heterogeneous
source amplitudes without FP32 norm underflow.

## Correctness and reproduction

The real-field C ABI tests cover both methods, full and odd solves,
reconstruction, source scales differing by 1e16, zero/exact/nonzero guesses,
collinear sources, a 3+2 partial tail, scalar agreement, per-column iteration
budgets and work accounting. A GCC-built client passes against the Apple
Clang dylib. The native Chroma regression performs 32 accepted checks,
including actual qprop and SuperbBLAS spin-3 colorvector sources, five RHS
submitted as device tiles of three plus two. All 19 standalone Objective-C++
programs pass syntax checks; C11 header and formatting checks pass. Fresh
Chroma configure accepts the batch API, rejects an old scalar-only library,
and accepts Metal-disabled builds.

After building/installing MetalClover and configuring Chroma with
`--with-metal=PREFIX`, build and run the native targets:

```sh
make -C mainprogs/tests t_metal_clover t_metal_clover_batch
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  mainprogs/tests/t_metal_clover gauge.lime action.xml 3,3,3,8
VECLIB_MAXIMUM_THREADS=1 OMP_NUM_THREADS=1 \
  mainprogs/tests/t_metal_clover_batch gauge.lime action.xml /tmp/batch12 24
```

The current command uses 24 color vectors (96 physical RHS), then 48 (192).
Use a fresh output prefix. The action XML follows the perambulator metadata
paths. The throughput harness measures vector qprop calls and is separate
from running the complete perambulator measurement.

[Structured measurements](results/batch-rhs12/summary.json),
[raw timings](results/batch-rhs12/run.csv), and
[validation records](results/batch-rhs12/validation/regression-summary.json)
record acceptance, widths, source hashes and work scope. Full logs and
resource samples are stored alongside them.
