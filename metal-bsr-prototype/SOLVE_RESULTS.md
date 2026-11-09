# Mixed-precision solve on M1 Max

Measured locally on November 3, 2026, on the same M1 Max (24 GPU cores, 32 GB, eight CPU worker tasks) as the multiply test. Connected to AC power, macOS 26.7.1. Each table entry is the median of three complete warm solves; backend order rotates between trials. Raw minimum/maximum variation can be recovered from the per-trial CSV.

All paths solve the same synthetic Hermitian positive-definite complex block-sparse system. The target is a **true FP64 relative residual ≤ 1e-10 for every RHS**, measured on the CPU over every component. Mixed paths use FP32 independent CG corrections with inner tolerance 0.002 and an FP64 accumulated solution. The all-FP64 CPU baseline uses CG at the final tolerance and verifies its true residual too. CPU kernels use multithreaded ARM NEON.

**Total solve time includes** input/residual conversion, FP32 iterations, all reductions and vector updates, Metal encoding, submission/waits, FP64 solution updates, and the full CPU FP64 residual checks. Allocation, operator upload, pipeline compilation, and manufactured-problem generation are excluded from this warm time and reported separately below. No kernel-only timing is used for the speedups.

## Complete solves: 4⁴ lattice

| RHS | CPU FP64 (ms) | CPU mixed (ms) | GPU mixed (ms) | GPU speedup vs faster CPU | GPU ms/RHS | GPU RHS/s |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1.11 | 2.81 | 9.98 | 0.11× | 9.978 | 100.2 |
| 12 | 4.41 | 9.04 | 13.37 | 0.33× | 1.114 | 897.9 |
| 48 | 11.99 | 18.81 | 28.02 | 0.43× | 0.584 | 1713.3 |
| 96 | 18.34 | 31.46 | 30.53 | 0.60× | 0.318 | 3144.3 |
| 192 | 32.12 | 66.17 | 40.72 | 0.79× | 0.212 | 4714.7 |
| 384 | 75.50 | 124.71 | 71.79 | 1.05× | 0.187 | 5349.2 |

## Complete solves: 8⁴ lattice

| RHS | CPU FP64 (ms) | CPU mixed (ms) | GPU mixed (ms) | GPU speedup vs faster CPU | GPU ms/RHS | GPU RHS/s |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 23.30 | 33.40 | 34.62 | 0.67× | 34.619 | 28.9 |
| 12 | 113.89 | 115.52 | 65.74 | 1.73× | 5.478 | 182.5 |
| 48 | 478.66 | 497.12 | 180.96 | 2.65× | 3.770 | 265.3 |
| 96 | 922.30 | 831.39 | 322.32 | 2.58× | 3.357 | 297.8 |
| 192 | 1877.31 | 1672.21 | 634.93 | 2.63× | 3.307 | 302.4 |
| 384 | 3836.78 | 3592.03 | 1287.30 | 2.79× | 3.352 | 298.3 |

## Convergence on 8⁴

Inner iterations below are the sum over refinement steps of the largest per-RHS iteration count, not a block-CG iteration count. GPU commands contain eight CG iterations each, so dispatched iterations also include padding after RHS have converged. The CPU stops immediately; the GPU's padded work is included in its timings.

| RHS | CPU FP64 iterations | CPU mixed corrections / inner iterations | GPU corrections / inner iterations / dispatched | Worst GPU true residual |
|---:|---:|---:|---:|---:|
| 1 | 36 | 4 / 55 | 4 / 55 / 80 | 6.01e-12 |
| 12 | 36 | 4 / 60 | 4 / 60 / 80 | 1.18e-11 |
| 48 | 36 | 4 / 60 | 4 / 60 / 80 | 1.16e-11 |
| 96 | 36 | 4 / 62 | 4 / 62 / 80 | 1.3e-11 |
| 192 | 36 | 4 / 62 | 4 / 62 / 80 | 1.46e-11 |
| 384 | 36 | 4 / 62 | 4 / 62 / 80 | 1.49e-11 |

Example convergence for 8⁴ / 48 RHS (first measured trial):

| Correction | CPU mixed max true residual | GPU mixed max true residual |
|---:|---:|---:|
| 1 | 0.00198 | 0.00198 |
| 2 | 3.76e-06 | 3.76e-06 |
| 3 | 6.68e-09 | 6.67e-09 |
| 4 | 1.18e-11 | 1.16e-11 |

## Where time goes and batching

GPU mixed time breakdown on the 8⁴ lattice (median of each component independently; rounded components need not sum to the median total):

| RHS | Inner solve (ms) | CPU FP64 residual (ms) | FP64 update (ms) | One-time GPU setup (ms) | Warm solve + setup (ms) |
|---:|---:|---:|---:|---:|---:|
| 1 | 32.47 | 1.64 | 0.36 | 1.14 | 35.76 |
| 12 | 54.04 | 10.04 | 0.84 | 2.59 | 68.33 |
| 48 | 131.22 | 44.03 | 2.36 | 32.99 | 213.94 |
| 96 | 233.52 | 78.22 | 4.26 | 20.25 | 342.57 |
| 192 | 450.30 | 165.09 | 8.65 | 56.43 | 691.36 |
| 384 | 907.79 | 348.99 | 17.43 | 245.42 | 1532.72 |

The highest observed GPU throughput on the larger lattice is at **192 RHS**, at 3.307 ms/RHS (302.4 RHS/s). This is the best of the tested batch sizes, not a proven optimum. Batch size should be judged by time per RHS rather than the wall time of an entire batch.

These are independent CG solves scheduled as a batch, with separate coefficients and convergence masks per RHS. There is no block Krylov method, shared Krylov basis, or multigrid preconditioner. More RHS can amortize dispatch overhead and improve utilization, but also increase memory traffic and the cost of CPU FP64 checks. Larger batches need not improve throughput indefinitely.

One-time setup includes workspaces, pipeline specialization and matrix/shared-buffer allocation and upload. A shared Metal library initialization of 37.83 ms is also excluded. Setup values come from the original allocations before warmup; adding them to warm medians is an accounting estimate, not a separately measured cold-start solve. CPU setup costs are in the raw CSV too.

## Validation and limits

Every measured solve passed full-component residual checks and a comparison against an independently generated known solution. The self-test covers nine small systems, zero RHS locking, odd RHS counts, and both SIMD remainder paths. A separate mass=0.01 run checks a higher condition-number bound. The original 30 multiply tests still pass after the utility-code refactor.

The operator is a periodic covariant Laplacian with fixed unitary color transports and spin permutations, plus mass=0.25. Its spectrum is bounded by [mass, mass+16]. It is a synthetic SPD problem with especially few distinct eigenvalues on very small lattices. It is not the Wilson/Dirac matrix in mgproton.cc. CG requires SPD; a real nonsymmetric operator needs a suitable solver and preconditioner. These results establish the tested mixed-precision mechanics and throughput, not production multigrid convergence or speedup. No Chroma/SuperbBLAS checkout was modified.

See [SOLVE.md](SOLVE.md) for algorithms and reproduction. Raw data:

- [results/solve.csv](results/solve.csv) and its `.trace.csv` file.
- [results/solve-extra.csv](results/solve-extra.csv) and its `.trace.csv` file.
