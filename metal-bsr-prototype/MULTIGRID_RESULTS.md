# Two-level clover multigrid on M1 Max

Measured locally on November 4, 2026, on the M1 Max (24 GPU cores, 32 GB, macOS 26.7.1, AC power). CPU runs use eight worker tasks and single-threaded Accelerate for the small coarse solve. Entries are medians of three complete solves from zero after a full warmup, with rotating method order. All methods share the same c_sw=1 clover operator and RHS within a case.

**The two-level implementation passed its numerical checks, but was slower than the existing BiCGStab baseline in these tests.** The coarse grid reduces FGMRES iteration counts; the tables below include the extra fine Schur applications and transfer work needed by each cycle, rather than judging performance by iterations alone.

## Complete solve times

BiCGStab is the existing even/odd clover baseline with FP32 elimination. FGMRES uses the same odd Schur operator and FP32 elimination; its preconditioner is either two clover-block smoothing steps or the two-level cycle. Full CPU FP64 residual checks and FP64 accumulation are retained for every method. Times include preparation, reconstruction, vector operations, reductions, small FGMRES QR solves, GPU submission/waits, and FP64 checks. Setup is excluded from warm solve times and reported separately.

| Mass | RHS | CPU BiCGStab (ms) | CPU two-level (ms) | Metal BiCGStab (ms) | Metal two-level (ms) | Metal MG gain vs BiCGStab | CPU MG / Metal MG |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 2265.67 | 3758.93 | 571.44 | 729.80 | 0.78× | 5.15× |
| 0.1 | 192 | 4442.49 | 7183.19 | 1323.80 | 1482.61 | 0.89× | 4.84× |
| 0.5 | 96 | 1681.54 | 2449.69 | 504.25 | 552.00 | 0.91× | 4.44× |
| 0.5 | 192 | 3643.07 | 5603.67 | 1073.36 | 1232.07 | 0.87× | 4.55× |

Gain is baseline time divided by multigrid time; below one means multigrid is slower. These comparisons use fresh timings from the same run, not historical clover results.

## Isolating the coarse-grid effect

| Mass | RHS | CPU smooth FGMRES (ms) | CPU two-level (ms) | Metal smooth FGMRES (ms) | Metal two-level (ms) | Metal coarse-grid gain |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 3160.50 | 3758.93 | 683.93 | 729.80 | 0.94× |
| 0.1 | 192 | 5907.98 | 7183.19 | 1373.18 | 1482.61 | 0.93× |
| 0.5 | 96 | 2020.39 | 2449.69 | 492.62 | 552.00 | 0.89× |
| 0.5 | 192 | 4197.65 | 5603.67 | 1046.74 | 1232.07 | 0.85× |

The smoother-only path has two damped block-Jacobi steps and one fine Schur application per preconditioner call. The two-level path inserts a coarse correction and recomputes the residual before post-smoothing, using two fine Schur applications per call. FGMRES then needs its own Schur application. A coarse-grid iteration reduction must offset that added cost.

## Convergence and actual operator work

Counts are `outer FP64 corrections / summed max per-RHS Krylov iterations / actual fine Schur calls / coarse calls`. BiCGStab and FGMRES iterations represent different work. Metal BiCGStab includes command-chunk padding; FGMRES has no such padding. FGMRES recomputes its true FP32 Schur residual at every restart, including when its small QR estimate suggests convergence.

| Mass | RHS | Method | Corrections / iterations / fine S / coarse | Worst full FP64 residual |
|---:|---:|---|---:|---:|
| 0.1 | 96 | cpu_bicgstab | 4 / 46 / 90 / 0 | 8.27e-12 |
| 0.1 | 96 | metal_bicgstab | 4 / 46 / 104 / 0 | 8.27e-12 |
| 0.1 | 96 | cpu_fgmres_smooth | 4 / 52 / 111 / 0 | 6.55e-12 |
| 0.1 | 96 | metal_fgmres_smooth | 4 / 52 / 111 / 0 | 6.55e-12 |
| 0.1 | 96 | cpu_fgmres_mg | 4 / 45 / 142 / 45 | 7.14e-12 |
| 0.1 | 96 | metal_fgmres_mg | 4 / 45 / 142 / 45 | 7.14e-12 |
| 0.1 | 192 | cpu_bicgstab | 4 / 47 / 90 / 0 | 8.7e-12 |
| 0.1 | 192 | metal_bicgstab | 4 / 47 / 112 / 0 | 8.71e-12 |
| 0.1 | 192 | cpu_fgmres_smooth | 4 / 53 / 113 / 0 | 6.46e-12 |
| 0.1 | 192 | metal_fgmres_smooth | 4 / 53 / 113 / 0 | 6.46e-12 |
| 0.1 | 192 | cpu_fgmres_mg | 4 / 45 / 142 / 45 | 6.64e-12 |
| 0.1 | 192 | metal_fgmres_mg | 4 / 45 / 142 / 45 | 6.64e-12 |
| 0.5 | 96 | cpu_bicgstab | 4 / 33 / 63 / 0 | 5.31e-11 |
| 0.5 | 96 | metal_bicgstab | 4 / 33 / 88 / 0 | 5.31e-11 |
| 0.5 | 96 | cpu_fgmres_smooth | 4 / 31 / 68 / 0 | 3.1e-12 |
| 0.5 | 96 | metal_fgmres_smooth | 4 / 31 / 68 / 0 | 3.1e-12 |
| 0.5 | 96 | cpu_fgmres_mg | 4 / 29 / 91 / 29 | 1.8e-12 |
| 0.5 | 96 | metal_fgmres_mg | 4 / 29 / 91 / 29 | 1.8e-12 |
| 0.5 | 192 | cpu_bicgstab | 4 / 33 / 63 / 0 | 6.7e-12 |
| 0.5 | 192 | metal_bicgstab | 4 / 33 / 88 / 0 | 6.7e-12 |
| 0.5 | 192 | cpu_fgmres_smooth | 4 / 31 / 68 / 0 | 3.1e-12 |
| 0.5 | 192 | metal_fgmres_smooth | 4 / 31 / 68 / 0 | 3.1e-12 |
| 0.5 | 192 | cpu_fgmres_mg | 4 / 29 / 91 / 29 | 1.68e-12 |
| 0.5 | 192 | metal_fgmres_mg | 4 / 29 / 91 / 29 | 1.68e-12 |

Full-system residual history for mass 0.1, 192 RHS:

| Correction | Metal BiCGStab | Metal smoother FGMRES | Metal two-level FGMRES |
|---:|---:|---:|---:|
| 1 | 0.0013 | 0.0013 | 0.0012 |
| 2 | 2.56e-06 | 2.39e-06 | 2.04e-06 |
| 3 | 4.53e-09 | 3.52e-09 | 3.84e-09 |
| 4 | 8.71e-12 | 6.46e-12 | 6.64e-12 |

## Hierarchy and amortization

The 8⁴ fine lattice has 24,576 complex odd-site components. 4⁴ aggregation leaves 16 aggregates with 12 coarse components each (six per chirality), giving a 192×192 coarse operator. Six seeded full-spin test vectors with constant and random parts are filtered by eight homogeneous block-Jacobi steps, split by chirality, and orthonormalized twice locally. Setup uses no benchmark RHS or known solution. P and restriction P† preserve gamma5. The FP64 Galerkin matrix is `A_c=P† S_o P`; its direct inverse is computed on the CPU once and rounded for FP32 CPU/GPU applications. One pre- and one post-smoothing step use omega=0.65. The same mass-specific hierarchy is reused at 96 and 192 RHS.

| Mass | RHS | Hierarchy setup (ms) | CPU MG warm ms/RHS | Metal MG warm ms/RHS | Metal MG RHS/s | Hierarchy break-even batches vs Metal BiCGStab |
|---:|---:|---:|---:|---:|---:|---|
| 0.1 | 96 | 664.07 | 39.156 | 7.602 | 131.5 | No warm-time gain |
| 0.1 | 192 | 664.07 | 37.412 | 7.722 | 129.5 | No warm-time gain |
| 0.5 | 96 | 467.38 | 25.518 | 5.750 | 173.9 | No warm-time gain |
| 0.5 | 192 | 467.38 | 29.186 | 6.417 | 155.8 | No warm-time gain |

Hierarchy setup includes filtering, local QR, probing P's columns, and dense coarse inversion. Problem generation and per-method workspace/pipeline setup are separate raw columns. BiCGStab and smoother-only solves do not need hierarchy setup; the field is repeated in their CSV rows only for case provenance. If the warm multigrid solve is slower, an unlimited supply of RHS cannot amortize away that per-solve overhead.

## Timing spread and FP64 cost

| Mass | RHS | CPU MG min–max (ms) | Metal MG min–max (ms) | Metal MG inner (ms) | Metal MG CPU FP64 checks (ms) |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 3572.47–3873.14 | 728.85–779.48 | 554.67 | 174.47 |
| 0.1 | 192 | 7016.70–8170.55 | 1466.93–1498.51 | 1114.72 | 349.40 |
| 0.5 | 96 | 2442.14–3077.00 | 536.60–598.21 | 374.67 | 168.60 |
| 0.5 | 192 | 4902.16–5629.98 | 1119.53–1640.62 | 811.47 | 370.09 |

Inner time includes FP32 RHS reduction/reconstruction. The Metal FGMRES vector, transfer, coarse-solve and fine-operator work stays on the GPU; only small Hessenberg coefficients/norms are read per Arnoldi step for CPU FP64 Givens rotations and triangular solves. Those synchronizations and calculations are inside the timing. The run used ordinary local system activity, without thermal or scheduler isolation.

## Validation and scope

All 72 measured solves met the full CPU FP64 true relative residual target of 1e-10 for every RHS, after even-site reconstruction. Maximum known-solution relative error was 3.06e-11. Small self-tests verify P†P, coarse gamma5 Hermiticity, the Galerkin identity, exact coarse correction on range(SP), CPU/Metal restriction and prolongation, V-cycle agreement, zero RHS locking, and all six complete solver paths. The actual 8⁴/4⁴-aggregate hierarchies also passed those operator identities at both benchmark masses, and the original 40-solve clover regression passed after the shared wrapper change.

This is a first two-level implementation on small positive-mass synthetic gauge fields. It has six local vectors per chirality and fixed smoothing/filtering parameters, not a tuned production adaptive near-null solve. Dense coarse inversion is intentionally limited to 512 coarse unknowns. Different aggregates, basis sizes, smoother/coarse accuracy and setup quality need separate measurements. The result does not establish production multigrid scaling or behavior near critical mass.

See [MULTIGRID.md](MULTIGRID.md) for equations, algorithms and reproduction.

Raw data: [results/multigrid.csv](results/multigrid.csv), [FP64 convergence traces](results/multigrid.csv.trace.csv).
