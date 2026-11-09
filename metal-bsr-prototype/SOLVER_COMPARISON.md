# Comparing clover Schur solvers on 16⁴

`metal_solver_compare` compares CPU and Metal versions of three mixed-precision
solvers on exactly the same synthetic clover-Wilson operator and RHS:

* BiCGStab on the odd clover Schur complement.
* FGMRES(8) with two clover-block Jacobi smoothing steps, weight 1.0.
* FGMRES(8) with the tuned sparse two-level cycle: twelve basis vectors per
  chirality, 4⁴ aggregates, one pre/post smoothing step at weight 1.0, and
  a sparse coarse BiCGStab correction with block-diagonal preconditioning,
  target preconditioned residual 0.1, and four-iteration cap.

All use c_sw=1, QDP gamma conventions, gauge strength 0.25, antiperiodic time,
positive bare masses 0.1/0.5, and the same odd Schur operator
`A_oo-H_oe*A_ee^-1*H_eo`. All inner arithmetic, reduction of the even/odd RHS,
and reconstruction are FP32 on the selected processor. The outer solution
accumulates in FP64. Every outer correction checks **all components of the full
original clover-Wilson residual in CPU FP64**, including the reconstructed even
sites. Every RHS must satisfy a relative residual of 1e-10. The inner target
is 0.002 for all methods. A zero RHS uses an absolute residual criterion.

The sparse hierarchy and its validation are described in
[SPARSE_MULTIGRID.md](SPARSE_MULTIGRID.md). This comparison changes neither the
operator nor the sparse multigrid settings. Smoother-only FGMRES uses the same
smoother and restart as the two-level method, with the coarse correction
omitted. BiCGStab uses the established unpreconditioned odd-Schur baseline.
Its GPU convergence check runs every iteration (`--bi-chunk 1`), avoiding
multi-iteration padding at this larger volume. Actual GPU Schur calls still
include the second half-step dispatch when RHS have already converged.

## Timing protocol and RHS pairing

96 RHS remain resident simultaneously. A logical 192-RHS batch consists of
two sequential 96-RHS tiles, sharing operator and solver state storage. Every
tile starts its solution from zero. CPU and Metal states are allocated one
method at a time so large Krylov stacks do not compete for the 32 GB Mac's
memory. The hierarchy is independent of all benchmark RHS and is shared across
methods/processors/trials for each mass. Baselines do not require hierarchy
setup; their hierarchy cost is recorded as zero.

Each trial solves global RHS 0–95 and then 96–191. The first tile supplies the
96-RHS timing; the sum of **both actually measured, distinct tiles** supplies
the 192-RHS timing. These are paired prefix measurements, rather than
independent repeated campaigns. No solve time is extrapolated. Hashes of the
fine component and global RHS index generate the manufactured solution, so
all methods see identical data regardless of execution order. The 96-RHS data
are exactly the prefix of the 192-RHS data.

Every method gets one complete resident-tile warmup, followed by three measured
trials. Method order rotates across mass/processor cases; processor order
reverses between masses. One method retains its workspace for its three
trials, then releases it before the next method is allocated. CPU uses eight
GCD worker tasks; Accelerate setup BLAS uses one internal thread. This runs
under ordinary local system activity, without scheduler or thermal isolation.

Complete solve time includes resets from zero, FP32 preparation/reconstruction,
fine/coarse operator work, smoothing, Krylov vectors and reductions, small CPU
QR, GPU submissions/waits, FP64 solution accumulation and full FP64 residual
checks. Known-solution error checking follows acceptance outside solve timing.
Manufactured RHS generation and applying M64 to generate b are outside solve
timing, with separate columns. Input construction is common to all methods.

Hierarchy, problem, inner workspace, and FP64 outer-workspace allocation are
measured separately. The same FP64 outer buffers are retained across all
methods, processors and trials for a mass/volume case. Their reset and input
movement are inside the solve timer. GPU allocation may reserve shared memory
lazily; cold first-use costs are covered by the warmup rather than claimed as
part of allocation timing. Setup tables are not complete cold-start totals.

## Operator work and convergence

One fine Schur call comprises two half-volume hops, the even clover inverse,
and the odd clover block. Counters include every dispatched fine application,
including preconditioner work and FGMRES residual refreshes. Iterations from
BiCGStab and FGMRES represent different work. Comparing iterations alone is
therefore insufficient.
Reported Krylov counts sum the maximum per-RHS iterations over outer
corrections. Coarse iteration counts sum the corresponding batch maxima over
preconditioner calls. These are batched iteration counts, not sums over all
96 individual RHS.

FGMRES checks its actual FP32 Schur residual on restart and retains the
preconditioned basis Z, permitting the inexact variable coarse correction.
The coarse solver's actual unpreconditioned FP32 residual is recomputed after
each correction and reported separately. It can exceed its preconditioned
target; reaching the iteration cap is allowed. Full FP64 convergence remains
mandatory for every method. CPU stops submitting coarse work on convergence;
Metal encodes the bounded coarse loop with per-RHS masks. Coarse GPU SpMV
counters include masked/padded dispatches. All of this work is timed.

CSV rows for 192 RHS sum work across both resident tiles; counts for a full
192-RHS batch must not be confused with a single 192-wide operator call.
Per-tile convergence ranges and full FP64 traces identify this explicitly.
Prefix rows share their first tile, so summing RHS counts over both logical
sizes double counts that tile. The report counts physical tile solves separately.

## Reproduction

```sh
make build/metal_solver_compare
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare --self-test
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \
  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass .1,.5 --csw 1 \
  --samples 3 --threads 8 --bi-chunk 1 --csv results/solver-compare16.csv
python3 solver_compare_summary.py
```

`make solver-compare-bench` runs the main comparison and generates its report.
`--backend cpu|metal|both` and `--methods bicgstab,smooth,sparse_mg` select a
subset for diagnostics. `--no-warmup` is intended for functional screening.
The small self-test uses two volumes, both masses, zero RHS, and both prefix
sizes, exercising all six paths. An allocation fix ensures Metal BiCGStab's
direction workspace remains full-sized after unused multiply-input arrays
are released; the established clover regression verifies compatibility.

Measurements, spread, convergence, setup amortization and raw data are in
[SOLVER_COMPARISON_RESULTS.md](SOLVER_COMPARISON_RESULTS.md). All code remains
in the standalone prototype, separate from Chroma and SuperbBLAS. These are
small positive-mass synthetic gauge problems; production ensemble and
critical-mass conclusions require separate validation.
The follow-up low-gap scan and negative bare-mass runs are documented in
[NEAR_CRITICAL.md](NEAR_CRITICAL.md). The driver now permits bare masses in
(-1,2], checks local clover factorization, and accepts `--coarse-iters`,
`--coarse-tol`, `--validate-only` and `--metal-first`. The original comparison
retains its positive masses and four-step coarse cap.
