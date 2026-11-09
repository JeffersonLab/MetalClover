# Two-level clover-Wilson preconditioner

This adds a two-level aggregation cycle and flexible GMRES to the standalone
prototype, preserving c_sw=1 clover blocks, Chroma/QDP conventions, odd-site
Schur elimination, FP32 CPU/Metal inner work, and **full CPU FP64 residual checks**.
See [MULTIGRID_RESULTS.md](MULTIGRID_RESULTS.md) for the original configuration
and [MULTIGRID_TUNING.md](MULTIGRID_TUNING.md) for the subsequent tuning and fresh
96/192-RHS comparisons. The executable defaults retain the original settings
for reproducibility; the tuned settings are explicit below.

## Fine operator and full stopping criterion

The fine multigrid operator is the retained odd clover Schur complement

```text
S_o = A_oo - H_oe A_ee^-1 H_eo,       H = -0.5 Dslash.
```

Each local A has two Hermitian 6×6 chiral blocks, with A_ee^-1 precomputed in
FP64 and rounded for FP32 use. The clover coefficient, gauge construction,
spin basis, and field-strength normalization are unchanged from [CLOVER.md](CLOVER.md).
No normal equations are formed.

All methods use the same mixed-precision outer correction loop. Every correction
reduces the current full residual to the odd system, approximately solves that
system, reconstructs the even correction, and accumulates into an FP64 solution.
The stopping criterion is **every RHS** satisfying

```text
||b - M64*x64||_2 / ||b||_2 <= 1e-10
```

on the full original clover-Wilson equation. The CPU evaluates every component
in FP64, including the clover term. Zero RHS use an absolute residual criterion.
Convergence is not accepted from an Arnoldi estimate or a reduced residual alone.

## Hierarchy setup

The default 8⁴ full lattice has 2,048 odd sites, or 24,576 complex fine components.
The setup aggregates full-space coordinate cubes of 4⁴ sites. Each contains 128
odd sites; there are 16 aggregates. Six basis vectors per chirality give 12 coarse
components per aggregate, for a **192×192** coarse operator. The tuned setting
uses twelve per chirality, giving a **384×384** coarse operator, while retaining
4⁴ aggregates. Basis size and all four aggregate dimensions are now runtime
parameters. Aggregate indexing uses full lattice coordinates, including time.

Setup is independent of the benchmark RHS and manufactured solution:

1. Generate six full-spin complex test vectors from constant random spin/color
   parts plus local random perturbations, using a fixed seed.
2. Apply eight homogeneous damped block-Jacobi filtering steps
   `v <- v - omega*A_oo^-1*S_o*v`, with omega=0.65, in FP64.
3. Split by the diagonal gamma5 eigenspaces and perform twice-reorthogonalized
   local Gram-Schmidt separately within each aggregate/chirality.
4. Use the resulting sparse prolongator P and its adjoint restriction R=P†.
   P has six nonzeros per fine spin/color component and obeys P†P=I and
   `gamma5_fine*P = P*gamma5_coarse`.
5. Form the FP64 Galerkin matrix `A_c=P†*S_o*P` by explicitly applying the fine
   operator to all prolongator columns. It retains coarse gamma5 Hermiticity.
6. Compute a pivoted dense inverse in FP64 on the CPU once, then round it and P
   for FP32 cycle application.

The same mass-specific hierarchy is reused across the 96- and 192-RHS cases.
Gauge links and clover coefficients are generated before RHS-dependent data,
so the hierarchy and runtime operators agree despite their different batch sizes.
Setup includes auxiliary matrix/vector construction and probing costs; it is
not claimed to be an optimized production setup.

This is smoothed-test-vector aggregation with fixed setup sweeps, not a
production adaptive near-null solve to an eigenvector tolerance. More local
vectors, other aggregates, and stronger setup may change the result. The dense
coarse implementation is limited to 512 unknowns. Larger coarse grids need
a sparse/iterative implementation rather than this direct inverse.

## Two-level cycle

For a fine residual r, one cycle returns an approximate correction z:

```text
z  = omega*A_oo^-1*r                          # pre-smooth
t  = r - S_o*z
e_c = A_c^-1 * P†*t                          # coarse correction
z += P*e_c
t  = r - S_o*z
z += omega*A_oo^-1*t                         # post-smooth
```

This is one pre-smoothing step, an exact-in-setup-precision coarse solve, and one
post-smoothing step. FP32 rounding makes the applied coarse solve approximate.
Every cycle needs two fine Schur applications. A smoother-only control omits
the coarse correction and the second residual computation, retaining two
block-Jacobi steps with one fine Schur application.

`--pre-smooth` and `--post-smooth` independently choose the number of steps.
For example, zero pre-smoothing starts with `z=P*A_c^-1*P†*r`; zero post-smoothing
returns immediately after coarse correction. With a coarse correction, a cycle
uses `pre+post` fine Schur applications. Without it, the same total number of
smoothing steps uses `pre+post-1` applications. FGMRES adds one Schur application
per Arnoldi step and one per residual refresh. Counters include all of these.
The tuned cycle retains one pre/post step and raises the **cycle** weight from
0.65 to 1.0. The independent **setup-filter** weight remains 0.65.

CPU cycle kernels use ARM NEON across RHS and eight GCD worker tasks. CPU coarse
application uses Accelerate complex GEMM with its internal thread count fixed
at one. Metal performs smoothing, restriction, dense coarse inverse application,
prolongation, and residual updates on the GPU, with barriers between dependent
dispatches. Local and coarse matrices are reused across RHS. The fine Schur
operator remains the existing generic hopping/local-block implementation.

## Flexible GMRES

The new right-preconditioned solver is restarted **FGMRES(8)**. It stores the
preconditioned basis Z explicitly, uses two-pass classical Gram-Schmidt, and
uses separate per-RHS complex Hessenberg systems and stopping masks. CPU and
GPU fine-vector arithmetic and reductions are FP32. Small Givens rotations and
triangular least-squares solves use the same CPU FP64 control code for both paths.

Metal keeps all fine Krylov vectors, preconditioned vectors, and cycle workspaces
resident. Each Arnoldi step returns only small reduced coefficients and norms
to the CPU; the full basis is not downloaded. Synchronization and small QR
control work are included in solve time. The code copies the current basis
vector into an operator buffer because the original Schur encoder expects
zero-offset buffers. These copies are timed and could be optimized later.

At every restart, including a predicted early stop, FGMRES recomputes the true
FP32 Schur residual `r=b32-S32*delta32`. The inner solve stops only when all RHS
meet the same 0.002 relative criterion used by the earlier mixed solver. If
needed, it restarts with that actual residual. The full outer loop still checks
the FP64 clover-Wilson equation. Krylov iteration limits, breakdown, nonfinite
values, failed true residuals, and known-solution mismatches cause errors.

FGMRES is used because it accommodates a general multigrid preconditioner for
this non-Hermitian operator and future variable coarse solves. The current
cycle itself is fixed. There is no sharing of a Krylov space across RHS:
these are independent solves executed in a batch.

## Six compared methods

The benchmark measures CPU and Metal versions of:

* Existing even/odd clover BiCGStab, without the new multigrid preconditioner.
* FGMRES with the smoother-only control.
* FGMRES with the two-level cycle.

All elimination/reconstruction is FP32 on the corresponding processor, while
FP64 accumulation and full residual checks remain on the CPU. All share the
same operator, RHS, outer tolerance, and inner tolerance. There is one full
warmup and three measured solves from zero per method, with rotating order.
The comparison does not rely on timings from earlier chats.

In the original configuration the coarse correction lowers FGMRES iteration counts,
but its extra Schur application and transfers outweigh that gain here. The
report therefore records **actual fine-operator calls** as well as iterations.
BiCGStab iterations and FGMRES iterations are not equivalent amounts of work.
Metal BiCGStab also includes its command-chunk padding; FGMRES has none.

Hierarchy setup is separated from warm solve timing, which includes everything
from zero/reset through final full FP64 residual acceptance. Problem generation,
workspace allocation, and compilation are also separate. Reusing the hierarchy
across many RHS amortizes setup; it cannot remove a slower warm cycle's per-RHS
cost. There is no favorable setup break-even against BiCGStab in the measured
configuration. The tuned report repeats this comparison rather than assuming
that reducing Krylov iterations or fine Schur calls guarantees a runtime gain.

## Build and reproduce

Run from this directory with Xcode tools and GPU access:

```sh
make
make mg-test
make mg-bench
python3 mg_summary.py results/multigrid.csv --output MULTIGRID_RESULTS.md
```

The benchmark defaults are:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/metal_multigrid \
  --sizes 8 --rhs 96,192 --mass .1,.5 --csw 1 --gauge-strength .25 \
  --aggregate 4 --setup-sweeps 8 --omega .65 --restart 8 \
  --samples 3 --threads 8 --csv results/multigrid.csv
```

A short complete solve comparison:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/metal_multigrid \
  --sizes 4 --rhs 1,3 --mass .25 --aggregate 2 --samples 1 \
  --csv results/mg-quick.csv
```

Validate the actual benchmark-size hierarchy without timing full solves:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/metal_multigrid \
  --sizes 8 --rhs 3 --mass .1,.5 --aggregate 4 --validate-only
```

Other options include `--tol`, `--inner-tol`, and `--periodic-time`. Aggregates
must be even and divide the even lattice extent. Restart is limited to 2–12;
the default memory guard limits a full vector to 12 million complex elements,
because CPU/GPU Krylov stacks retain many half-volume vectors. The comparison
holds workspaces for several methods and is not a minimum-memory benchmark.

The tuned configuration is:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/metal_multigrid \
  --sizes 8 --rhs 96,192 --mass .1,.5 --csw 1 --gauge-strength .25 \
  --aggregate-shape 4,4,4,4 --vectors 12 --setup-sweeps 8 --setup-omega .65 \
  --pre-smooth 1 --post-smooth 1 --omega 1 --restart 8 \
  --samples 3 --threads 8 --csv results/multigrid-tuned.csv
```

`make mg-tune-screen` reproduces three sequential screening stages.
`make mg-tuned-bench` runs numerical checks and the fresh original/tuned
comparison, then regenerates `MULTIGRID_TUNING.md`. These commands require GPU
access. The latter alternates configuration order across mass/RHS cases, warms
each method, and records three measured solves from zero. BiCGStab is measured
with the original configuration; its algorithm does not depend on the hierarchy.
The tuned smoother-only control uses the same stronger smoothing weight as MG.

Additional options: `--vectors` (1–24 per chirality), `--aggregate-shape`
(four even dimensions), `--setup-sweeps` (0–96), `--setup-omega` (0–1),
`--omega` (0–1.5), `--pre-smooth`, `--post-smooth` (1–6 total), and `--methods`
(comma-separated driver method names). `--no-warmup` is for convergence
screening, not reported warm timings. Dense coarse size remains limited to 512.

## Validation and scope

Self-tests verify local P†P, coarse gamma5 Hermiticity, the Galerkin identity
`A_c*y=P†*S_o*P*y`, and inverse consistency. They also verify that coarse
correction solves a residual in range(S_o P), compare CPU and Metal transfers
and cycles, exercise zero RHS locking, and run all six complete solvers with
full FP64 residual/known-solution checks. `--validate-only` checks the same
hierarchy identities at the actual 8⁴ benchmark geometry.

The final checks passed for both benchmark-size hierarchies and the earlier
40-solve clover self-test. Logs are saved as
`results/multigrid-hierarchy-validation.log` and `results/clover-mg-regression.log`.

The known solution is used only for post-solve validation, not to create P.
Every full residual component and every RHS is checked. The matrices remain
synthetic, positive bare mass, and small volume. This establishes a working
two-level preconditioner and its measured costs, not production multigrid
scaling or behavior near a critical mass. The follow-up sweep tunes aggregation,
basis size, filtering, smoothing count/weight, and restart on these fields.
It does not establish optimal parameters for a physical ensemble or a larger
lattice, and the dense coarse solver remains a scalability limit.

Files: `multigrid.hpp` (hierarchy and cycles), `fgmres.hpp` (CPU/Metal solver),
`multigrid.metal` (device transfers, coarse solve, Krylov operations),
`mg_checks.hpp` (independent checks), `multigrid.mm` (comparison driver), and
`mg_summary.py` (report). All preexisting standalone drivers remain available.
