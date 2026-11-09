# M1 Max Metal block-sparse benchmark

Measured locally on November 3, 2026: Apple M1 Max, 24-core GPU, 10-core CPU (8 performance + 2 efficiency), 32 GB unified memory; macOS 26.7.1 (25G241), connected to AC power. CPU parallel runs request eight worker tasks; Accelerate internal threads are limited to one.

Source data: [results/benchmark.csv](results/benchmark.csv). Times below are median **wall milliseconds per multiply**, including command encoding/submission/wait for Metal. Initial allocation, data copies, and pipeline compilation are excluded. All data remains resident.

CPU best means the fastest measured CPU implementation in that case, chosen from single-worker Accelerate, eight-worker Accelerate, and eight-worker NEON. GPU batch amortizes command-buffer overhead across independent repeated multiplies; GPU single submits and waits for one multiply. Speedup uses CPU best divided by the corresponding GPU time.

## Dense 12×12 complex blocks

| Lattice | RHS | CPU best (ms) | CPU method | GPU single (ms) | GPU batch (ms) | Single speedup | Batch speedup |
|---|---:|---:|---|---:|---:|---:|---:|
| 4⁴ | 1 | 0.091 | NEON 8 | 0.286 | 0.045 | 0.32× | 2.01× |
| 4⁴ | 4 | 0.116 | BLAS 8 | 0.444 | 0.040 | 0.26× | 2.93× |
| 4⁴ | 12 | 0.368 | BLAS 8 | 0.276 | 0.076 | 1.33× | 4.86× |
| 4⁴ | 48 | 1.091 | NEON 8 | 0.563 | 0.186 | 1.94× | 5.86× |
| 8⁴ | 1 | 0.733 | BLAS 8 | 0.546 | 0.121 | 1.34× | 6.04× |
| 8⁴ | 4 | 0.803 | BLAS 8 | 0.327 | 0.140 | 2.46× | 5.72× |
| 8⁴ | 12 | 2.267 | BLAS 8 | 1.116 | 0.288 | 2.03× | 7.87× |
| 8⁴ | 48 | 17.483 | NEON 8 | 2.011 | 1.667 | 8.69× | 10.49× |
| 12⁴ | 1 | 3.651 | BLAS 8 | 1.806 | 0.691 | 2.02× | 5.28× |
| 12⁴ | 4 | 4.013 | BLAS 8 | 2.002 | 0.711 | 2.00× | 5.65× |
| 12⁴ | 12 | 11.439 | BLAS 8 | 1.594 | 1.383 | 7.18× | 8.27× |
| 12⁴ | 48 | 79.693 | NEON 8 | 8.553 | 8.446 | 9.32× | 9.44× |
| 16⁴ | 1 | 10.499 | NEON 8 | 2.490 | 2.292 | 4.22× | 4.58× |
| 16⁴ | 4 | 11.880 | BLAS 8 | 2.808 | 2.503 | 4.23× | 4.75× |
| 16⁴ | 12 | 35.092 | BLAS 8 | 4.982 | 4.843 | 7.04× | 7.25× |
| 16⁴ | 48 | 249.435 | NEON 8 | 28.234 | 28.014 | 8.83× | 8.90× |

## 3×3 color blocks with four-spin permutation

| Lattice | RHS | CPU best (ms) | CPU method | GPU single (ms) | GPU batch (ms) | Single speedup | Batch speedup |
|---|---:|---:|---|---:|---:|---:|---:|
| 4⁴ | 1 | 0.060 | NEON 8 | 0.250 | 0.022 | 0.24× | 2.79× |
| 4⁴ | 4 | 0.073 | NEON 8 | 0.253 | 0.022 | 0.29× | 3.32× |
| 4⁴ | 12 | 0.109 | BLAS 8 | 0.287 | 0.041 | 0.38× | 2.68× |
| 4⁴ | 48 | 0.261 | BLAS 8 | 0.244 | 0.044 | 1.07× | 5.88× |
| 8⁴ | 1 | 0.242 | NEON 8 | 0.297 | 0.074 | 0.81× | 3.28× |
| 8⁴ | 4 | 0.368 | NEON 8 | 0.284 | 0.048 | 1.30× | 7.65× |
| 8⁴ | 12 | 1.000 | NEON 8 | 0.676 | 0.144 | 1.48× | 6.93× |
| 8⁴ | 48 | 3.280 | BLAS 8 | 0.786 | 0.624 | 4.17× | 5.26× |
| 12⁴ | 1 | 1.003 | NEON 8 | 0.498 | 0.096 | 2.01× | 10.46× |
| 12⁴ | 4 | 1.670 | NEON 8 | 0.411 | 0.206 | 4.07× | 8.10× |
| 12⁴ | 12 | 5.055 | NEON 8 | 1.932 | 0.697 | 2.62× | 7.26× |
| 12⁴ | 48 | 19.010 | BLAS 8 | 3.237 | 3.235 | 5.87× | 5.88× |
| 16⁴ | 1 | 3.476 | NEON 8 | 1.073 | 0.373 | 3.24× | 9.33× |
| 16⁴ | 4 | 5.342 | NEON 8 | 1.264 | 0.678 | 4.23× | 7.87× |
| 16⁴ | 12 | 15.417 | NEON 8 | 3.004 | 2.260 | 5.13× | 6.82× |
| 16⁴ | 48 | 53.583 | BLAS 8 | 10.188 | 10.191 | 5.26× | 5.26× |

## Interpretation and limits

Across these 32 synthetic cases, amortized GPU speedups range from 2.01× to 10.49× versus the fastest measured CPU method. Single-submit GPU execution is faster in 26/32 cases. These are case-specific measurements, not a prediction of full-solver speedup.

All numerical checks passed. The largest observed absolute error was 9.04e-07; the largest relative L2 error was 2.67e-07 against an independent double-precision reference. Small benchmark cases are fully checked; large cases use seeded samples and boundary sites. The separate 30-case self-test fully checks every output, including irregular sparsity and empty rows. An additional [full-output check](results/full-check.log) verifies all 589,824 outputs for each operator on an 8⁴ lattice with 12 RHS. That separate check is not included in the primary timing table.

The largest per-implementation matrix/input/output working set was 1233.3 MiB. The process holds separate CPU and GPU copies.

This demonstrates a standalone FP32 complex multiply using the relevant sparse structures on this M1 Max. It does not benchmark SuperbBLAS itself or implement a complete backend. Actual multigrid performance also depends on reductions, tensor rearrangements, precision conversion, CPU/GPU dependencies, and numerical convergence. Batched throughput is especially optimistic for workflows that must return to the CPU between every operation.

The primary run used 5 samples per backend and case. Timing spread and batch sizes are preserved in the CSV; GPU latency for small cases is sensitive to OS scheduling and submission overhead. This was a local run with ordinary system activity, not an isolated laboratory measurement.

See [README.md](README.md) for the equations, timing definitions, validation coverage, and reproduction commands.
