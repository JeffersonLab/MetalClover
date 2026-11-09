# Mixed-precision solve extension

This adds a complete small solve to the standalone prototype. Everything remains
in this directory; Chroma and SuperbBLAS are neither linked nor modified. See
[SOLVE_RESULTS.md](SOLVE_RESULTS.md) for measured timings, convergence, and RHS
throughput on the M1 Max. The older multiply benchmark remains available.

## Build and reproduce

From this directory:

```sh
make
make solve-test
make solve-bench
python3 solve_summary.py results/solve.csv --output SOLVE_RESULTS.md
```

Default benchmark: `4⁴` and `8⁴` lattices; 1, 12, 48, 96, 192, and 384 RHS;
three measured solves per backend, plus one full warmup. For a quick example:

```sh
./build/metal_solve --sizes 4 --rhs 12,48 --samples 1 \
  --csv results/solve-quick.csv
```

The saved main study was run in two parts:

```sh
./build/metal_solve --rhs 1,12,48,192 --csv results/solve.csv
./build/metal_solve --rhs 96,384 --csv results/solve-extra.csv
python3 solve_summary.py results/solve.csv results/solve-extra.csv \
  --output SOLVE_RESULTS.md
```

Additional saved checks:

```sh
./build/metal_solve --self-test --csv results/solve-self-test.csv
./build/metal_solve --sizes 4 --rhs 7 --mass 0.01 --samples 1 \
  --csv results/solve-stress.csv
```

Other options: `--threads 8`, `--chunk 8`, `--mass 0.25`, `--tol 1e-10`, and
`--inner-tol 0.002`. Shader paths are relative to the prototype directory.
CPU numerical kernels use ARM NEON and GCD, not Accelerate BLAS, in this solve.
The solve therefore does not require a BLAS threading environment variable.
The older multiply benchmark still requires `VECLIB_MAXIMUM_THREADS=1`.
GPU access must be available; a restrictive automation sandbox can hide devices.

## Test problem and precision

The block structure is the same compressed 3-color / 4-spin structure used in the
multiply prototype. There are nine blocks per lattice site, with periodic shifts
in four directions. The matrix is a covariant Laplacian plus a positive mass:

```text
A = mass * I + sum_mu (2 I - T_mu - T_mu^H)
```

`T_mu` is a periodic forward shift combined with a fixed dense unitary color
matrix, a cyclic spin permutation, and a complex unit phase. The backward edge
uses their adjoints. This guarantees that A is Hermitian positive definite for
positive mass, with spectrum in `[mass, mass+16]`. The default mass is 0.25,
giving a condition-number upper bound of 65. The mass=0.01 stress check has
upper bound 1601. Small periodic lattices have few distinct eigenvalues and can
converge unusually quickly; these are not realistic gauge configurations.

The color matrix for direction mu is `F diag(1, i^(mu+1), (-1)^(mu+1)) F^H`,
where F is the unitary three-point Fourier matrix and mu ranges from zero to
three. The spin shift is mu+1 and its unit phase is i^mu. Diagonal site blocks
are `(mass+8) I`. Both directions and all blocks are explicitly stored; the
CPU/GPU kernels do not exploit the operator's translational invariance.

The authoritative matrix coefficients, manufactured solution, RHS, accumulated
solution, and true residuals are FP64. The GPU uses separately rounded FP32
matrix coefficients; its residual checks are **not** performed on an FP32-rounded
matrix promoted back to FP64. The RHS is manufactured as `b = A64 * x_known`
using a seeded double-precision solution. Each benchmark RHS is nonzero and
different. Self-tests also include an exactly zero RHS.

This problem is deliberately SPD so conjugate gradient is valid. The actual
Wilson/Dirac application may require GCR, GMRES, BiCGStab, or another suitable
method, with a preconditioner. This extension does not establish convergence for
that operator and does not implement multigrid.

## Solver paths

Three implementations solve the same system, starting at zero:

1. **cpu_fp64**: independent-per-RHS CG entirely on the CPU in FP64, targeting
   1e-10, followed by a full true-residual check. If that check misses the target,
   the same outer correction mechanism can refine it further.
2. **cpu_mixed**: FP32 CG corrections on the CPU, double-precision accumulation
   and true-residual checks. It uses the same correction tolerance and stopping
   policy as the GPU mixed path.
3. **metal_mixed**: FP32 CG corrections on Metal, with FP64 accumulation and
   true-residual checks on the CPU.

Each mixed solve performs iterative refinement (defect correction):

```text
x64 = 0; r64 = b64
repeat:
    approximately solve A32 * delta32 = float(r64) with CG
        each RHS stops when its recursive residual is <= 0.002 of its initial norm
    x64 += double(delta32)
    r64 = b64 - A64*x64                    # full CPU computation
    stop only when every RHS has ||r64||2 / ||b64||2 <= 1e-10
```

A RHS that meets the true tolerance is locked by supplying a zero correction
RHS in subsequent outer steps. A zero original RHS uses an absolute residual
criterion and stays zero. CG has a 256-iteration inner limit and refinement has
a 12-step limit; breakdown, nonfinite values, or failure to converge cause a
nonzero exit, not a timing presented as success.

CPU FP32 dot products use FP32 accumulation, matching the GPU's arithmetic
precision. CPU FP64 dot products and all true-residual norms use FP64. Reduction
orders differ between CPU and GPU, so convergence need not be bit-for-bit
identical. Inner correction limits are intentionally looser than final accuracy.

## GPU execution and batching

Matrix, correction, recursive residual, search direction, matrix-product output,
and reduction scratch remain in shared Metal buffers during each inner solve.
There are GPU kernels for the sparse matrix product, dot-product partial sums,
reduction/alpha, correction and residual update, reduction/beta, and search
direction update. No recursive residual or CG scalar is computed on the CPU.

The GPU executes `--chunk` CG iterations per command buffer, with barriers between
dependent dispatches, then waits so the CPU can inspect the per-RHS active flags.
The default chunk is eight. Once an individual RHS converges, its direction is
zeroed and updates are masked, but unmasked sparse products and reductions can
still perform padded work until the command chunk ends. Those costs are included
in total runtime and reported as dispatched iteration counts. The CPU CG stops
as soon as all RHS converge, without this padding.

This is **batched independent CG**, not block CG. Each RHS has its own rho, alpha,
beta, iteration count, and convergence mask. Adding RHS amortizes command
submission and scalar-control overhead and provides more parallel work. It
does not reduce the mathematical work per RHS through a shared Krylov space.

With an unlimited RHS supply, choose batch size by throughput (`RHS / total
time`) rather than the time for a whole batch. On this M1 Max and the tested 8⁴
problem, throughput is roughly flat at 96–384 RHS, with the best measured median
at 192. The 192 versus 384 difference is small relative to ordinary run-to-run
variation; there is no evidence here that a larger batch is inherently better.
For a larger lattice, retune the batch size under the memory limit. The current
CLI caps a case at 24 million complex vector elements to keep the numerous
FP32/FP64 work vectors manageable on this 32 GB machine.

## Timing, validation, and files

The reported **complete solve time** starts with resetting the zero solution and
copying the RHS into the residual. It includes casts to FP32, CG iterations and
all vector/reduction work, Metal encoding/submission/waits, CPU FP64 correction
updates, and all full FP64 true-residual checks. Initial problem generation,
manufacturing b, computing fixed b norms, matrix upload, workspace allocation,
pipeline compilation, and library initialization are outside that warm time.
Setup costs are separately recorded; workspaces and uploaded operators are reused
across the repeated solves. Each timed solve restarts from zero; it is not a warm
start from an earlier solution. The same RHS data is reused across timing trials.

Known-solution error is computed after the timed solve, as a separate validation.
The true-residual check is inside the timing and covers **every component of
every RHS**, with no sampling. The known-solution error is required to be finite
and within a condition-number-scaled tolerance. All solves use a manufactured
solution, not convergence of a recursive FP32 residual alone, to check accuracy.
Self-tests additionally verify reverse-edge adjoints, spin inverses, color
unitarity, and the positive-definiteness energy bound independently of CG.

The CSV contains one row per measured solve: total, inner, residual, and update
times; setup costs; outer/inner/dispatched iterations; worst true residual;
known-solution relative error; milliseconds per RHS; and RHS per second. A
companion `.trace.csv` records true-residual convergence after each correction.
`solve_summary.py` creates the report using only Python's standard library.

`solve.mm` holds the CG benchmark driver; `solve_common.hpp` holds the operator,
CPU SIMD products, CPU CG, Metal CG scheduling, and outer refinement shared with
the Wilson extension. `solve.metal` adds the CG kernels to the
existing `kernels.metal` sparse multiply at runtime. Shared helpers formerly in
`main.mm` now live in `bsr_common.hpp`; the original multiply benchmark and its
numerical checks are preserved.

The later [Wilson / BiCGStab extension](WILSON.md) replaces the SPD test operator
and inner solver while retaining the FP64 outer refinement/checking machinery.
