# Clover odd-Schur / Metal FP32 elimination results

Measured locally on November 4, 2026, on the M1 Max (24 GPU cores, 32 GB, eight CPU worker tasks), macOS 26.7.1, connected to AC power. Three complete solves from zero are measured per method after one full warmup, with rotating method order. All methods in a case share identical gauge links, clover blocks, and RHS.

The operator uses **c_sw=1.0**, QDP's gamma matrices with diagonal gamma5, and **two Hermitian 6×6 chiral blocks at every site**. The retained odd system is exactly `S_o = A_oo - H_oe A_ee^-1 H_eo`, without scalar mass division or Schur normalization. H includes the Wilson `-1/2 Dslash` factor. The links are seeded synthetic SU(3) matrices at strength 0.25, with periodic space and antiperiodic fermionic time boundaries.

Chroma's normalization is explicit: `F_AH=(Q-Q†)/8` from four oriented leaves and `C=-(c_sw/2) sum_(mu<nu) gamma_mu gamma_nu F_AH`. This agrees with the inspected Chroma packed-component formulas. The convention, local inverses, and full clover solves were validated before implementing the Metal elimination change.

Every solve is accepted only when **every RHS has full CPU FP64 true relative residual ≤ 1e-10 after even-site reconstruction**. FP64 accumulated solutions and checks use the original FP64 clover-Wilson coefficients. Mixed inner BiCGStab tolerance is 0.002.

## Complete solve runtime

All times include input conversion, elimination, Schur BiCGStab, reductions and vector updates, GPU submission/waits, reconstruction, FP64 accumulation, and full FP64 residual checks. Allocation, clover construction/inversion, and compilation are outside warm time and recorded separately. CPU best is the faster of the three measured CPU methods, not a tuned production Dslash baseline.

| Mass | Lattice | RHS | CPU best (ms) | CPU mixed FP32 elim. (ms) | Metal CPU-FP64 elim. (ms) | Metal device-FP32 elim. (ms) | Move gain | CPU best / Metal new | New ms/RHS |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 8⁴ | 96 | 3335.00 | 3335.00 | 910.07 | 674.18 | 1.35× | 4.95× | 7.023 |
| 0.1 | 8⁴ | 192 | 6800.84 | 7089.20 | 2054.47 | 1565.66 | 1.31× | 4.34× | 8.154 |
| 0.5 | 8⁴ | 96 | 2362.22 | 2362.22 | 852.50 | 575.38 | 1.48× | 4.11× | 5.994 |
| 0.5 | 8⁴ | 192 | 4704.36 | 4704.36 | 1718.19 | 1266.87 | 1.36× | 3.71× | 6.598 |

## All methods and convergence

Names ending in `host64` use CPU FP64 RHS reduction and reconstruction around FP32 inner solves. `cpu_mixed_fp32` performs all elimination math in CPU FP32. `metal_mixed_fp32` does that math, packing/scattering, and inner initialization on Metal. All paths retain CPU FP64 true-residual checks. Counts are `corrections / summed maximum inner iterations / dispatched iterations`; Metal pads to four-iteration chunks.

| Mass | RHS | Method | Total (ms) | Corrections / inner / dispatched | Worst true residual | Known-solution error |
|---:|---:|---|---:|---:|---:|---:|
| 0.1 | 96 | cpu_fp64 | 5142.58 | 1 / 41 / 41 | 6.59e-11 | 3.9e-10 |
| 0.1 | 96 | cpu_mixed_host64 | 3338.64 | 4 / 46 / 46 | 8.28e-12 | 2.93e-11 |
| 0.1 | 96 | cpu_mixed_fp32 | 3335.00 | 4 / 46 / 46 | 8.27e-12 | 2.93e-11 |
| 0.1 | 96 | metal_mixed_host64 | 910.07 | 4 / 46 / 52 | 8.27e-12 | 2.95e-11 |
| 0.1 | 96 | metal_mixed_fp32 | 674.18 | 4 / 46 / 52 | 8.27e-12 | 2.95e-11 |
| 0.1 | 192 | cpu_fp64 | 11232.22 | 1 / 41 / 41 | 6.57e-11 | 4.01e-10 |
| 0.1 | 192 | cpu_mixed_host64 | 6800.84 | 4 / 47 / 47 | 8.7e-12 | 3.05e-11 |
| 0.1 | 192 | cpu_mixed_fp32 | 7089.20 | 4 / 47 / 47 | 8.7e-12 | 3.05e-11 |
| 0.1 | 192 | metal_mixed_host64 | 2054.47 | 4 / 47 / 56 | 8.64e-12 | 3.05e-11 |
| 0.1 | 192 | metal_mixed_fp32 | 1565.66 | 4 / 47 / 56 | 8.71e-12 | 3.06e-11 |
| 0.5 | 96 | cpu_fp64 | 2441.12 | 1 / 19 / 19 | 6.59e-11 | 2.12e-10 |
| 0.5 | 96 | cpu_mixed_host64 | 2376.26 | 4 / 33 / 33 | 5.31e-11 | 1.52e-11 |
| 0.5 | 96 | cpu_mixed_fp32 | 2362.22 | 4 / 33 / 33 | 5.31e-11 | 1.52e-11 |
| 0.5 | 96 | metal_mixed_host64 | 852.50 | 4 / 33 / 44 | 5.31e-11 | 1.52e-11 |
| 0.5 | 96 | metal_mixed_fp32 | 575.38 | 4 / 33 / 44 | 5.31e-11 | 1.52e-11 |
| 0.5 | 192 | cpu_fp64 | 5167.51 | 1 / 19 / 19 | 6.63e-11 | 2.21e-10 |
| 0.5 | 192 | cpu_mixed_host64 | 4976.38 | 4 / 33 / 33 | 6.7e-12 | 9.7e-12 |
| 0.5 | 192 | cpu_mixed_fp32 | 4704.36 | 4 / 33 / 33 | 6.7e-12 | 9.69e-12 |
| 0.5 | 192 | metal_mixed_host64 | 1718.19 | 4 / 33 / 44 | 6.7e-12 | 9.69e-12 |
| 0.5 | 192 | metal_mixed_fp32 | 1266.87 | 4 / 33 / 44 | 6.7e-12 | 9.69e-12 |

## Time breakdown of the Metal change

Component medians are independent and may not sum to total. The inner column subtracts preparation/reconstruction subfields from raw `inner_ms`, avoiding double counting.

| Mass | RHS | Elimination | Schur inner (ms) | Preparation (ms) | Reconstruction (ms) | CPU FP64 residual (ms) | FP64 update (ms) |
|---:|---:|---|---:|---:|---:|---:|---:|
| 0.1 | 96 | CPU FP64 | 444.03 | 121.88 | 121.12 | 214.08 | 6.20 |
| 0.1 | 96 | Metal FP32 | 369.40 | 17.43 | 17.76 | 205.82 | 5.93 |
| 0.1 | 192 | CPU FP64 | 871.43 | 258.93 | 331.27 | 588.72 | 17.96 |
| 0.1 | 192 | Metal FP32 | 782.28 | 122.95 | 61.37 | 571.55 | 17.10 |
| 0.5 | 96 | CPU FP64 | 366.07 | 122.53 | 116.70 | 237.50 | 6.46 |
| 0.5 | 96 | Metal FP32 | 315.41 | 27.62 | 9.17 | 211.79 | 5.77 |
| 0.5 | 192 | CPU FP64 | 684.14 | 231.93 | 253.96 | 505.53 | 18.60 |
| 0.5 | 192 | Metal FP32 | 609.59 | 136.62 | 61.06 | 418.59 | 10.85 |

Metal preparation still includes a CPU FP64-to-FP32 cast of the full residual into shared memory, because Metal has no native FP64 arithmetic. The reduced odd RHS stays on the GPU; it is not downloaded to initialize BiCGStab. The FP32 full correction is read only after GPU reconstruction completes, then accumulated in CPU FP64. End-to-end timing includes these boundaries.

## Full FP64 residual history: mass 0.1, 192 RHS

| Correction | Metal + CPU FP64 elimination | Metal FP32 elimination |
|---:|---:|---:|
| 1 | 0.0013 | 0.0013 |
| 2 | 2.56e-06 | 2.56e-06 |
| 3 | 4.53e-09 | 4.53e-09 |
| 4 | 8.64e-12 | 8.71e-12 |

## Batching and timing spread

| Mass | RHS | New Metal min–max (ms) | New Metal RHS/s | Workspace setup (ms) | Problem setup (ms) |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 604.05–677.09 | 142.4 | 125.43 | 265.94 |
| 0.1 | 192 | 1425.88–1611.06 | 122.6 | 138.55 | 373.63 |
| 0.5 | 96 | 557.85–609.40 | 166.8 | 12.85 | 247.24 |
| 0.5 | 192 | 1160.83–1335.47 | 151.6 | 135.30 | 449.74 |

Metal library initialization was 48.44 ms. Problem setup includes links, clover construction, FP64 inverses, and manufactured RHS generation; workspace setup includes per-method buffers and pipelines. Shared resources are charged to each method separately, so setup columns should not be summed across methods.

Compare time per RHS and the timing spread when choosing a batch. These are independent batched BiCGStab solves, not a shared or block Krylov solver. The run used ordinary system activity and was not thermally or scheduler isolated.

## Validation and limits

All 60 measured solves passed full FP64 residual and known-solution checks. Stage 1 validated Chroma normalization and packed formulas, Hermitian chiral blocks, independent local inverses, gauge covariance, direct full clover and odd Schur products, and complete solves before device elimination changed. The final self-test adds direct device RHS/reconstruction checks and all five solver paths on eight small cases, including odd RHS, a zero RHS, free field, and c_sw=0. An analytic quantized Abelian field independently fixes the plaquette normalization and orientation.

The FP64 inverses are precomputed once and rounded for FP32 application; there is no GPU factorization claim. Clover/local matrices are reused across RHS. Field strengths are built once from gauge links; adding clover does not automatically reduce Dslash link traffic. The data remain synthetic, isotropic, and positive bare mass with no multigrid, production ensemble, or negative-critical-mass validation. Generic kernels are used on both CPU and GPU. The historical scalar-diagonal Wilson results are not an apples-to-apples baseline for the new clover operator.

See [CLOVER.md](CLOVER.md) for exact equations and reproduction.

Raw data: [results/clover.csv](results/clover.csv), [convergence traces](results/clover.csv.trace.csv).
