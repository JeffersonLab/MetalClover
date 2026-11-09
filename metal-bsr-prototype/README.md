# Standalone Metal complex block-sparse prototype

This directory is independent of Chroma and SuperbBLAS. It neither includes nor
links SuperbBLAS and does not modify that checkout. It benchmarks **single
precision complex** `y = A x` on synthetic periodic four-dimensional lattices.

The measured results and interpretation are in [RESULTS.md](RESULTS.md).

The real-field performance comparison with **Chroma's original reliable
mixed-precision clover BiCGStab inverter**, GPU BiCGStab, and reused GPU
multigrid is documented in [PERFORMANCE_REAL_FIELDS.md](PERFORMANCE_REAL_FIELDS.md).

The **mixed-precision solve extension**, CPU comparisons, and large-RHS batching
study are described in [SOLVE.md](SOLVE.md), with measured results in
[SOLVE_RESULTS.md](SOLVE_RESULTS.md).

The subsequent **Wilson / complex BiCGStab** extension is documented in
[WILSON.md](WILSON.md), with 96–192 RHS results in
[WILSON_RESULTS.md](WILSON_RESULTS.md).

Even/odd Schur-complement preconditioning and its complete CPU/Metal comparison
are in [EVENODD.md](EVENODD.md) and [EVENODD_RESULTS.md](EVENODD_RESULTS.md).

The **clover extension** uses Chroma/QDP conventions, two 6×6 chiral blocks per
site, the odd-site Schur complement, and Metal FP32 RHS reduction/reconstruction.
See [CLOVER.md](CLOVER.md) and [CLOVER_RESULTS.md](CLOVER_RESULTS.md).

The two-level aggregation preconditioner and flexible-GMRES comparison are in
[MULTIGRID.md](MULTIGRID.md) and [MULTIGRID_RESULTS.md](MULTIGRID_RESULTS.md).
The subsequent basis/aggregation/smoother tuning and fresh 96–192 RHS timings
are in [MULTIGRID_TUNING.md](MULTIGRID_TUNING.md).

The sparse iterative coarse-grid extension and 16⁴ comparison are in
[SPARSE_MULTIGRID.md](SPARSE_MULTIGRID.md) and
[SPARSE_MULTIGRID_RESULTS.md](SPARSE_MULTIGRID_RESULTS.md). It supports 96
resident RHS and bounded-memory processing of larger RHS batches.

The direct 16⁴ comparison against BiCGStab and smoother-only FGMRES is in
[SOLVER_COMPARISON.md](SOLVER_COMPARISON.md), with fresh measurements in
[SOLVER_COMPARISON_RESULTS.md](SOLVER_COMPARISON_RESULTS.md).

The negative bare-mass spectral scan and near-critical-like solver comparison
are in [NEAR_CRITICAL.md](NEAR_CRITICAL.md) and
[NEAR_CRITICAL_RESULTS.md](NEAR_CRITICAL_RESULTS.md).

The compact 32⁴ operator, third multigrid level, and RHS tiling are described in
[THREE_LEVEL.md](THREE_LEVEL.md). The full comparison against BiCGStab and
smoother-only FGMRES is reported in [THREE_LEVEL_RESULTS.md](THREE_LEVEL_RESULTS.md).

[CHROMA_ADAPTER.md](CHROMA_ADAPTER.md) describes the standalone NERSC/Chroma
bridge, the FP64 operator match, and input and output for batches of spinors.
[CHROMA_ADAPTER_RESULTS.md](CHROMA_ADAPTER_RESULTS.md) gives the real-field
CPU and Metal measurements.

## Real anisotropic clover field

`Makefile.production` builds the rectangular-lattice driver `production_solve`
and the native Chroma tools `production_export` and `production_verify`.
The exporter reads the valence action and times from the perambulator metadata,
checks LIME/QIO input, applies Chroma's FP64 gauge reunitarization, and performs
the spatial stout steps in FP64. Hopping anisotropy and the spatial/electric
clover coefficients are applied separately. The Metal driver uses FP32 odd
Schur solves in a refinement loop with full CPU FP64 residuals, one resident RHS, and one
three-level hierarchy reused across the four source spins.

Build with a Chroma source snapshot matching the installed library; an unrelated
development checkout can have incompatible headers. On this machine:

```sh
make -f Makefile.production production \
  CHROMA_CONFIG=/Users/admin/Downloads/macos-chroma/install/chroma/bin/chroma-config \
  CHROMA_SOURCE=/Users/admin/Downloads/qcd-macos-share-20261005/chroma.devel
```

Native parameter parsing and stout smearing are compiled locally in FP64;
the Chroma and SuperbBLAS source and installation trees are read only.
`production_export GAUGE_LIME PERAM_XML EIGEN_MOD PREFIX` writes the physical
stout links, anisotropic parameters, a full Chroma operator probe, and the stored
v0 vectors in the interchange defined by `production_bridge.hpp`. The solver
requires no Chroma calls after import:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/production_solve build/production24 results/production24 \
  > results/production24-solve.log 2>&1
OMP_NUM_THREADS=8 ./build/production_verify build/production24 results/production24 \
  > results/production24-chroma-verify.log 2>&1
python3 compare_production.py results/production24 \
  --chroma-residuals results/production24.chroma-residuals.csv \
  --solve-log results/production24-solve.log --action results/anisotropic24-action.json
```

Use a new output prefix for a fresh run. `--resume` rechecks saved full solutions
against the imported operator before reusing them. `--initial-prefix PREFIX`
loads available solutions as initial guesses and corrects their full residuals.
The measured 24³×256 run used the raw-link diagnostic as an initial guess for
spin 0; spins 1–3 started from zero. Exact commands, hashes and solver settings
are in `results/production24-provenance.json`. Timings in the comparison JSON
separate the hierarchy and solves; they exclude import and checkpoint I/O.

The comparison covers all 768 `(eigen_sink,eigen_source)=(0,0)` complex elements
at source time 176 and sink times 176–223. The historical perambulator has an
input target of 1e-8 with `SilentFail=true`, so its achieved residual is unknown.
The supplied FP32 gauge file's historical processing is also unconfirmed.
Report its measured projection discrepancy separately from the new solutions'
full FP64 residuals; the standalone driver enforces 1e-10 on every accepted solve.

The completed run passed the independent native Chroma checks for all four
sources, with a worst full relative residual of **1.752e-11**. Across the 768
perambulator elements the relative 2-norm discrepancy is **3.438e-5**, and the
maximum absolute discrepancy is **1.012e-5**. Its cause remains unresolved.
See [the complete comparison](results/production24.comparison.json),
[all element differences](results/production24.comparison.csv), and
[the native full residuals](results/production24.chroma-residuals.csv).

## Fresh 12³×96 correctness test

The test in `results/correctness12` uses the new LIME configuration and its own
header action: mass -0.0743, xi0=4.3, nu=1.265,
cR=1.58954576550983, cT=0.902866347885228, and two spatial stout steps with
rho=0.14. All eight fermion-action entries in that header agree. The fermion
boundary is `1 1 1 -1`; `_per` denotes antiperiodic time in these ensemble filenames.
No eigenvectors or historical perambulator values are used.

A dense, unit-norm complex Z4 source with seed 20261107 is generated independently
of the operator and stored once. Chroma's original `UnprecCloverLinOp` and
`InvBiCGStab`, compiled in FP64, give a fresh reference solution. Metal BiCGStab
and three-level MG-preconditioned FGMRES both start from zero on the same stored
RHS, use FP32 asymmetric odd Schur solves, accumulate corrections in FP64, and
require the **full** FP64 relative residual to be at most 1e-12. Native Chroma
then rereads the original LIME field and independently checks every saved solution.

| Path | Independent full Chroma FP64 relative residual |
| --- | ---: |
| Chroma FP64 BiCGStab | 4.6772e-14 |
| Metal BiCGStab + FP64 refinement | 8.4464e-15 |
| Metal MG-FGMRES + FP64 refinement | 1.8346e-14 |

The Metal BiCGStab and MG solution vectors differ by **4.4484e-13 relative**,
comparing all 1,990,656 complex components. Relative differences from the Chroma
reference are 2.0061e-12 and 1.7163e-12. Residual tolerance and forward solution
error are distinct quantities. The independent native stout-state and full
operator probes match the export exactly; the standalone FP64 full operator
matches at 1.6828e-16, and CPU/Metal FP32 Schur products at 7.5364e-8.

The following are the measured workflow commands; the exporter and Metal driver
require unused output paths when starting a new run. The Python preparation and
report scripts describe this particular recorded 12³×96 test.

```sh
python3 correctness12_prepare.py
make -f Makefile.production build/production_export build/correctness_chroma build/correctness_compare \
  CHROMA_CONFIG=/Users/admin/Downloads/macos-chroma/install/chroma/bin/chroma-config \
  CHROMA_SOURCE=/Users/admin/Downloads/qcd-macos-share-20261005/chroma.devel
gauge12=$(python3 -c 'import json; print(json.load(open("results/correctness12/input.json"))["gauge"])')
OMP_NUM_THREADS=8 ./build/production_export "$gauge12" results/correctness12/action.xml - build/correctness12 \
  > results/correctness12/export.log 2>&1
OMP_NUM_THREADS=8 VECLIB_MAXIMUM_THREADS=1 ./build/correctness_chroma \
  "$gauge12" results/correctness12/action.xml build/correctness12 results/correctness12/test reference \
  > results/correctness12/chroma-reference.log 2>&1
VECLIB_MAXIMUM_THREADS=1 ./build/correctness_compare build/correctness12 results/correctness12/test \
  > results/correctness12/metal.log 2>&1
OMP_NUM_THREADS=8 VECLIB_MAXIMUM_THREADS=1 ./build/correctness_chroma \
  "$gauge12" results/correctness12/action.xml build/correctness12 results/correctness12/test verify \
  > results/correctness12/chroma-verify.log 2>&1
python3 correctness12_report.py
```

[summary.json](results/correctness12/summary.json) contains the full residual
and solution comparisons, action parameters, iteration
traces, precision choices, timings and file hashes. These are correctness results;
the original Chroma reference uses a full unpreconditioned system, while the Metal
paths use even/odd preconditioning, so their iteration counts are not equivalent.

## 16³×128 correctness repeat

The same harness also accepts an explicit gauge file, output directory, and
header mass. The 16³×128 configuration contains both light and strange actions;
this repeat selects the light mass **-0.0840**, supported by six identical
header-action entries. Its cR=1.58932722549812 and cT=0.902783591772098,
xi0=4.3, nu=1.265, two rho=0.14 spatial stout steps, and `1 1 1 -1` fermion BC
are preserved from that header. The source seed remains 20261107, with one
newly sized, unit-normalized stored RHS shared by all solver paths.

```sh
python3 correctness12_prepare.py \
  --gauge /Users/admin/Downloads/download/prop_tests/szscl21_16_128_b1p50_t_x4p300_um0p0840_sm0p0743_n1p265_per_cfg_1000.lime \
  --mass=-0.0840 --out-dir results/correctness16
# Export, native reference, and native verification use the same commands as
# above with build/correctness16 and results/correctness16/test as prefixes.
VECLIB_MAXIMUM_THREADS=1 ./build/correctness_compare \
  build/correctness16 results/correctness16/test 4,4,4,4 --defer-reference \
  > results/correctness16/metal.log 2>&1
python3 correctness12_report.py \
  --folder results/correctness16 --bridge-prefix build/correctness16
```

The optional last Metal-driver argument selects the fine aggregate dimensions;
4×4×4×4 divides this geometry and leaves even first-coarse extents for the
second 2×2×2×2 aggregation. Solver precisions, the three-level construction,
and the full 1e-12 FP64 acceptance criterion are unchanged. Exact command arrays
and environment settings are recorded in
[workflow.json](results/correctness16/workflow.json).
The measured repeat overlaps the CPU reference and Metal solves. The optional
`--defer-reference` skips the initial load of the reference solution while it is
being computed; it changes no solver mathematics or residual threshold. Native
Chroma verification and the final complete-vector comparison still run after
both processes finish, and the report requires all three native residuals.

The completed light-quark repeat passed all independent full Chroma FP64 checks:

| Path | Full Chroma FP64 relative residual |
| --- | ---: |
| Chroma FP64 BiCGStab | 4.6718e-14 |
| Metal BiCGStab + FP64 refinement | 1.2750e-14 |
| Metal MG-FGMRES + FP64 refinement | 2.5298e-14 |

Across all **6,291,456 complex components**, the Metal BiCGStab and MG solution
vectors differ by **1.7082e-12 relative**. Their differences from the independent
Chroma reference are 2.5154e-12 and 3.0775e-12. The 1e-12 requirement applies to
the full residual; relative solution error is a separate quantity. The
standalone FP64 operator matches Chroma at 1.9859e-16, and CPU/Metal FP32 Schur
products match at 7.5875e-8. Native stout-state and full-operator probes match
the export exactly. See [the complete report](results/correctness16/summary.json)
and [independent native checks](results/correctness16/test.chroma-checks.csv).

## Build and run

Requires an Apple Silicon Mac, macOS 13.3 or newer, and Xcode command-line build
tools. The shader is compiled at runtime through Metal, so the optional offline
Metal compiler is not required. No Python packages or other libraries are needed.
Run these commands from this directory:

```sh
make
make test
make bench
python3 summarize.py results/benchmark.csv --output RESULTS.md
```

The executable must have access to the GPU. A restricted automation sandbox may
hide Metal devices; if that happens, run it in a normal Terminal session.

A shorter, fully checked run:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/metal_bsr \
  --sizes 4,8 --rhs 1,12 --samples 5 --full-check \
  --csv results/quick.csv
```

Other options: `--threads 8`, `--target-ms 20`, `--mode dense12`,
`--mode color3_spin4`, `--shader kernels.metal`. `--sizes` means the extent **L** of
an `L × L × L × L` lattice, not the number of sites. The default run covers four
extents (4, 8, 12, 16) and four RHS counts (1, 4, 12, 48), for both operators.

## Operators and layout

Each regular row has nine blocks: the site itself and its two periodic neighbors
in each of four directions. Row pointers and column indices are explicit CSR
arrays; the kernel does not infer neighbors from lattice coordinates.

* **dense12**: arbitrary dense 12×12 complex blocks.
  `y[i,c,r] = sum(e in row i, k) A[e,c,k] * x[col[e],k,r]`.
* **color3_spin4**: a dense 3×3 color block multiplied by a four-spin permutation
  with phases in `{1, -1, i, -i}`. Only the 3×3 block is stored, not a dense 12×12
  expansion. This models the structural specialization found in SuperbBLAS's
  `bsr_kron_3x3_4x4perm` code.
  `y[i,c,s,r] = sum(e,k) A[e,c,k] * phase[e,s] * x[col[e],k,perm[e,s],r]`.

Complex numbers are interleaved real/imaginary FP32 pairs; RHS is the fastest
dimension, then spin (where present), color, and lattice site. Matrix blocks use
row-major storage. Both operators have 12 complex output components per site per
RHS. Values and phases are reproducible from a fixed seed; these are **not** gauge
configurations, a Wilson operator, or a complete multigrid solve. Self-tests also
use irregular sparsity and entirely empty rows.

## Implementations

* `accelerate_1`: single worker, Accelerate complex GEMV/GEMM per sparse block.
* `accelerate_parallel`: the same operation split across eight GCD worker tasks.
  This follows the per-block BLAS strategy in SuperbBLAS's non-MKL CPU path;
  it is a separate implementation, **not a timing of SuperbBLAS itself**. In spin
  mode, the input spin permutation/phase is materialized in scratch before GEMM.
* `neon_parallel`: eight GCD tasks, ARM NEON over groups of four RHS values, scalar
  remainder, register accumulation across blocks. One RHS uses the scalar path.
* `metal_latency`: one multiply per command buffer, including CPU encoding,
  submission, GPU execution, and the completion wait.
* `metal_batch`: repeated independent multiplies on resident data, encoded in one
  compute encoder/command buffer, with a buffer memory barrier between dispatches.
  Includes CPU encoding, submission, GPU execution, and one completion wait,
  divided by the number of multiplies. This measures amortized throughput,
  **not the latency of a solver iteration**. It reuses the same input each time.

Accelerate's internal thread count is held at one with
`VECLIB_MAXIMUM_THREADS=1`, avoiding nested threading. Eight worker tasks is a
requested partition count, not CPU-core affinity or a guarantee of eight
simultaneously running performance cores. The fastest of all three CPU timings
is used for comparisons. Apple may schedule work differently during other runs.

The GPU assigns one complex output to each thread; there are no atomics. The
shader specializes RHS, color, and spin dimensions with function constants.
It uses safe Metal math, FP32 FMA, shared buffers, and 128-thread groups. It does
not use matrix instructions, FP16, or Apple-specific assembly. The C++ compiler
targets Apple M1 with `-O3 -ffp-contract=fast`, without `-ffast-math`.

## Timing and correctness

Data generation, shader/pipeline compilation, buffer allocation, and initial
copies are outside the multiply timings. CSV columns `upload_ms` (allocation plus
initial copies) and `pipeline_ms` record those costs separately. No input/output
transfers are included per multiply: GPU buffers are CPU-visible shared memory
and remain resident. CPU implementations likewise reuse preallocated input,
output, and scratch buffers; spin packing and zeroing needed on each apply are
included in CPU timing.

Each backend is warmed up. CPU and batched GPU repeat counts are calibrated to
approximately 20 ms per sample, capped at 128 repetitions; slow cases use one
repeat. Five samples are taken by default. The CSV records their median, minimum,
maximum, and repeats. GPU latency uses one dispatch per sample. CPU and GPU order
alternate by case. Device timestamps are also recorded (`device_ms`) but are
**not** substituted for wall time in the comparisons.

Useful GFLOP/s counts eight real operations per complex color multiply-add.
Permutation/phase handling, data movement, and packing are excluded from that
FLOP count. Spin mode counts work in the compressed representation, not a dense
12×12 expansion. `working_set_mib` includes one matrix, index/phase arrays, and one
input/output pair; the program holds separate CPU and GPU copies and additional
scratch, so actual process memory is larger.

An independent scalar **double-precision** implementation provides the numerical
oracle. All three CPU implementations and the GPU are checked. Every output is
checked in the 30 self-tests; these include RHS 1, 3, 4, 7, and 12, periodic
duplicate neighbors, missing blocks, and empty rows. Benchmark cases check every
output when there are at most 65,536 outputs; larger cases check 4,096 seeded
random positions plus all components/RHS at four selected sites. The CSV reports
the actual number checked. `--full-check` checks every output for every case.

Checks require finite results, relative L2 error <= 5e-6 and maximum absolute
error <= 2e-5 against the oracle. The batched GPU output is checked again after
timing. These checks concern this multiply only; they do not establish solver
convergence or double-precision suitability.

The saved `results/full-check.csv` and `.log` additionally verify every output of
both 8⁴-lattice / 12-RHS operators (589,824 complex outputs each). This is a
separate validation run; the primary timing table uses `results/benchmark.csv`.

## Files

* `main.mm`: multiply benchmark driver.
* `bsr_common.hpp`: shared deterministic problems, CPU implementations, reference
  checks, Metal setup, and multiply timing helpers.
* `kernels.metal`: generic and spin-structured complex block-sparse shader.
* `solve.mm`, `solve.metal`: mixed-precision solve extension; see `SOLVE.md`.
* `solve_common.hpp`: shared CPU SIMD products, solver types, and FP64 refinement.
* `wilson.mm`, `bicgstab.hpp`, `bicgstab.metal`: Wilson operator and BiCGStab.
* `wilson_common.hpp`: shared Wilson construction and independent stencil checks.
* `evenodd.mm`, `evenodd.hpp`, `evenodd.metal`: packed even/odd Schur solver.
* `clover.mm`, `clover_common.hpp`, `clover_checks.hpp`, `clover.metal`: clover
  construction, local block kernels, odd Schur solver, and device elimination.
* `multigrid.mm`, `multigrid.hpp`, `fgmres.hpp`, `mg_checks.hpp`,
  `multigrid.metal`: two-level hierarchy, cycles, CPU/Metal FGMRES, and checks.
* `mg_tune_sweep.py`, `mg_tuned_bench.py`, `mg_tuned_summary.py`: reproducible
  multigrid tuning, fresh comparisons, and the [tuning report](MULTIGRID_TUNING.md).
* `sparse_multigrid.mm`, `sparse_coarse.hpp`, `sparse_coarse.metal`,
  `sparse_checks.hpp`: sparse Galerkin assembly, iterative coarse CPU/Metal
  solves, large-volume benchmark, and validation.
* `run_sparse_study.py`, `sparse_summary.py`: sequential regression/benchmark
  workflow and the 16⁴ results report.
* `solver_compare.mm`, `sparse_rhs.hpp`, `solver_compare_summary.py`: shared
  deterministic RHS and the paired 96/192-RHS comparison of all three solvers.
* `critical_probe.mm`, `nearcritical_summary.py`: full FP64 spectral diagnostics
  and the low-gap negative-mass comparison report.
* `compact_clover.hpp`, `compact_clover.metal`, `compact_solve.hpp`,
  `three_level.hpp`, `three_level.mm`: compact stencil, memory reuse and the
  recursive three-level CPU/Metal solver at 32⁴.
* `compact_probe.mm`, `threelevel_summary.py`: fixed-memory spectral screening
  and audited totals for the complete 96/192-RHS sweeps.
* `results/`: raw CSV and text logs from the local runs.
* `summarize.py`: regenerates the results report using Python's standard library.

The prototype is intentionally limited to forward `y = A x`, FP32 complex, one
GPU, and a fixed operator/input during each timing sample. It does not implement
conjugate transpose, general alpha/beta, MPI, precision conversion, tensor
permutations, solver dependencies, or SuperbBLAS's allocation/stream interfaces.
Those remain work for a full backend.
