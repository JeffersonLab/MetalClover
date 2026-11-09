# Sparse iterative coarse grid for clover-Wilson

`metal_sparse_multigrid` replaces the global dense coarse inverse in the
standalone solver with an explicit block-sparse Galerkin matrix and an inexact,
block-preconditioned complex BiCGStab solve. The old small-lattice dense driver
remains available as an independent reference. Chroma and SuperbBLAS checkouts
are unchanged.

The operator keeps c_sw=1, the QDP basis with diagonal gamma5, two local 6×6
clover blocks per fine site, and

```
S_o = A_oo - H_oe A_ee^-1 H_eo.
```

Fine FGMRES and all CPU/Metal inner work use FP32. The solution accumulates in
CPU FP64. After every outer correction, even sites are reconstructed and **all
components of the original full clover-Wilson residual are evaluated in FP64**.
Every RHS must satisfy `||b-M64*x64||/||b|| <= 1e-10`. A zero RHS uses an absolute
criterion. A coarse residual or a GMRES estimate cannot accept the full solve.

## Sparse setup

The basis uses the previous tuned settings: twelve vectors per chirality,
4⁴ aggregates, eight FP64 test-vector filtering sweeps at weight 0.65, and
local twice-reorthogonalized QR. The cycle uses one pre/post clover-block
smoothing step at weight 1.0 and FGMRES restart 8.

On 16⁴ there are 32,768 odd sites, 393,216 complex fine components, and 256
aggregates. Each aggregate has 24 coarse components: **6,144 coarse unknowns**.
The Schur operator's two-hop paths generate a coarse pattern with 33 neighboring
blocks per aggregate, including the diagonal, for **8,448 blocks of size 24×24**.
No coupling is discarded using a numerical threshold.

Setup never forms a global dense coarse matrix or all prolongator columns:

1. At each even site, form the small blocks `T=H_eo*P`, grouped by the source
   aggregate. Only aggregates touching that even site are stored.
2. Apply the local FP64 inverse `A_ee^-1` and gamma5 to those blocks.
3. Assemble, by disjoint coarse rows,
   `A_c = P† A_oo P - gamma5_c T† gamma5_e A_ee^-1 T`.
   This uses `P† H_oe = gamma5_c (H_eo P)† gamma5_e`.
4. Store the result as block CSR in FP64 and FP32. Invert each *local* 24×24
   diagonal block in FP64 for block-Jacobi preconditioning. There is no global
   coarse inverse or global factorization.

At 16⁴ the FP32 sparse matrix occupies 37.125 MiB. Its FP64 copy occupies
74.25 MiB; the basis, local inverses, and sparse indices are additional.
Independent dense column probing on small lattices verifies the assembly,
including entries outside the sparse pattern. At 16⁴, random full FP64
`A_c*y` and `P†*S_o*P*y` checks verify the actual benchmark geometry.

## Iterative coarse correction

Each coarse correction starts from zero and solves the left-preconditioned
non-Hermitian system

```
D_c^-1 A_c e = D_c^-1 r_c,
```

where `D_c` contains the local 24×24 blocks. BiCGStab is suitable here; no normal
equations or positive-definiteness assumption is introduced. The target is a
0.1 relative **preconditioned** residual, with an explicit iteration cap. The
actual unpreconditioned FP32 residual `r_c-A_c*e` is recomputed after every
coarse solve and reported separately. It need not be below 0.1, and reaching
the cap is allowed: this is an inexact preconditioner for flexible GMRES.
Breakdowns and nonfinite residuals cause errors; only the full FP64 stopping
criterion accepts the final solution.

CPU coarse sparse products use ARM NEON across RHS and eight worker tasks.
Metal performs sparse products, local preconditioning, BiCGStab scalar
reductions, per-RHS stopping masks, and updates entirely on the GPU, inside the
fine-cycle command encoder. It encodes up to the configured cap; masked RHS
skip further sparse arithmetic. CPU observes small diagnostics when the
existing fine Arnoldi step completes. These extra dispatches, reductions,
synchronization, and diagnostics are included in solve time. Raw GPU coarse
SpMV counters include masked/padded dispatches, whereas CPU stops submitting
products after convergence.

FGMRES retains its preconditioned vectors because the capped/tolerance-based
coarse solve gives a variable preconditioner. Its actual FP32 Schur residual
is recomputed on restart. The inner fine target remains 0.002.

## Large RHS batches and memory

The sparse coarse matrix removes the dense coarse setup/solve bottleneck;
fine Krylov storage still scales with lattice volume and resident RHS count.
CPU and Metal workspaces are allocated **sequentially**, so they do not retain
two large Krylov stacks simultaneously. Unused random multiply-input arrays
are released. The driver exposes `--tile-rhs` to bound memory.

Logical RHS batches are processed as whole resident tiles, reusing the hierarchy,
operator, pipelines, and inner workspace. `--tile-rhs` must divide the requested RHS
count. All RHS are distinct, generated from component/global-RHS indices;
96 RHS are a prefix of the corresponding 192-RHS problem, independent of tile
width. CPU and Metal receive the same operator and RHS. Each tile starts from
zero and checks every RHS. Batch residual/error fields are maxima over tiles.

Reported batch solve time sums complete tile solve times, including FP64
checks, RHS reduction/reconstruction, reset, and GPU overhead. Synthetic RHS
construction and its application to generate b are outside the solve timer
and have their own column. Hierarchy, problem, and workspace setup are also
separate. The workspace column covers inner-solver/elimination allocation;
FP64 outer buffers are created once per logical case outside the solve timer,
and their allocation time was not collected. Reset/copy work inside solve is
timed. Thus the setup columns do not give a complete cold-start total.
Trace files record global RHS offsets, so tiling is explicit.
A timing at 192 logical RHS must not be described as 192 simultaneously resident
RHS when the reported tile width is smaller.

## Build, validate, reproduce

Run from this directory with Xcode command-line tools and GPU access:

```sh
make build/metal_sparse_multigrid
VECLIB_MAXIMUM_THREADS=1 ./build/metal_sparse_multigrid --self-test
VECLIB_MAXIMUM_THREADS=1 ./build/metal_sparse_multigrid \
  --sizes 16 --rhs 3 --tile-rhs 3 --validate-only
```

The full comparison accepts `--rhs 96,192`, `--tile-rhs`, `--mass .1,.5`,
`--coarse-iters`, `--coarse-tol`, `--samples 3`, and `--backend cpu|metal|both`.
`--csw 1` is accepted and other values are rejected for this comparison.
`--smoother-only` runs the same fine FGMRES cycle without a coarse correction.
`--no-warmup` is intended for screening, not final warm timings.

Self-tests compare sparse assembly with the old dense Galerkin oracle, check
coarse gamma5 Hermiticity, compare CPU/Metal sparse products and strict coarse
solves, exercise zero RHS, and solve complete clover-Wilson systems on CPU and
Metal. The 16⁴ geometry also has independent Galerkin/coarse-solve validation.
The existing dense multigrid and clover drivers provide regression checks for
the shared generic BiCGStab and coarse-callback interfaces.

This remains a single-device prototype on synthetic gauge fields and positive
bare masses. The sparse implementation removes the 512-unknown dense limit;
production adaptive setup, distributed sparse coarse operators, deeper levels,
and critical-mass behavior require further work. Measurements and exact
commands are in [SPARSE_MULTIGRID_RESULTS.md](SPARSE_MULTIGRID_RESULTS.md).
The follow-up direct comparison with BiCGStab and smoother-only FGMRES uses
[SOLVER_COMPARISON.md](SOLVER_COMPARISON.md) and
[SOLVER_COMPARISON_RESULTS.md](SOLVER_COMPARISON_RESULTS.md).
