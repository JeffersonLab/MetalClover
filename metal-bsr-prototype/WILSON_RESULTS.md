# Wilson / BiCGStab results on M1 Max

Measured locally on November 3, 2026: M1 Max, 24 GPU cores, 32 GB memory, macOS 26.7.1, AC power. CPU runs request eight GCD workers and use ARM NEON. Each entry is the median of three complete solves from a zero initial guess, after a full warmup. Backend order rotates between trials.

This is the **unpreconditioned Wilson operator**, not the earlier SPD Laplacian: `(m+4)I - 0.5 Dslash`, with spin projectors, site-dependent SU(3) links, periodic space and antiperiodic time. The links are seeded synthetic near-identity matrices, not a QCD gauge ensemble. The benchmark uses gauge strength 0.25 and positive bare masses 0.1 and 0.5.

All three paths use **complex BiCGStab directly on D**, with independent recurrence coefficients and stopping masks for every RHS. Mixed paths use FP32 corrections (inner tolerance 0.002), FP64 accumulation, and full CPU FP64 true-residual checks. The CPU FP64 baseline uses the final tolerance directly. Every RHS must have true relative residual ≤ 1e-10. No normal equations or CG are used in these Wilson solves.

## Complete solve runtime

Times include all BiCGStab vector operations and reductions, operator applications, Metal encoding/submission/waits, precision conversion, FP64 solution updates, and full CPU FP64 residual checks. Reused buffer allocation, matrix upload, and compilation are excluded from warm solve time and itemized separately.

| Mass | Lattice | RHS | CPU FP64 (ms) | CPU mixed (ms) | Metal mixed (ms) | Speedup vs faster CPU | Metal ms/RHS | RHS/s |
|---:|---|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 8⁴ | 96 | 6693.56 | 4287.48 | 889.93 | 4.82× | 9.270 | 107.9 |
| 0.1 | 8⁴ | 192 | 13046.18 | 8027.42 | 1757.28 | 4.57× | 9.153 | 109.3 |
| 0.5 | 8⁴ | 96 | 3422.78 | 2964.90 | 651.11 | 4.55× | 6.782 | 147.4 |
| 0.5 | 8⁴ | 192 | 7179.16 | 6026.93 | 1415.50 | 4.26× | 7.372 | 135.6 |

## Convergence

Iteration counts sum the maximum per-RHS BiCGStab count at each correction step. Each iteration needs up to two Wilson applications. The GPU executes commands in groups of four iterations and includes padded work after RHS converge; the CPU can stop immediately, including after the first half of an iteration. Padded work is included in the Metal timings.

| Mass | RHS | CPU FP64 corrections / iterations | CPU mixed corrections / iterations | Metal corrections / iterations / dispatched | Worst Metal true residual |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 1 / 72 | 4 / 86 | 4 / 85 / 92 | 1.31e-11 |
| 0.1 | 192 | 1 / 71 | 4 / 87 | 4 / 87 / 92 | 1.34e-11 |
| 0.5 | 96 | 1 / 37 | 4 / 59 | 4 / 58 / 64 | 1.26e-11 |
| 0.5 | 192 | 1 / 38 | 4 / 61 | 4 / 61 / 72 | 1.32e-11 |

True-residual history: mass 0.1, 8⁴, 192 RHS, first measured trial.

| Correction | CPU mixed | Metal mixed |
|---:|---:|---:|
| 1 | 0.002 | 0.002 |
| 2 | 3.98e-06 | 3.98e-06 |
| 3 | 7.44e-09 | 7.46e-09 |
| 4 | 1.44e-11 | 1.34e-11 |

## Runtime components and variation

Component medians are independent and need not sum to the median total. Setup is the original pre-warmup allocation and upload cost. It is not a measured cold-start total.

| Mass | RHS | CPU mixed min–max (ms) | Metal total min–max (ms) | FP32 inner (ms) | CPU FP64 checks (ms) | FP64 update (ms) | Metal setup (ms) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 4285.81–5043.05 | 878.06–941.52 | 710.78 | 164.60 | 4.61 | 20.79 |
| 0.1 | 192 | 7940.93–8195.67 | 1686.75–1832.28 | 1398.06 | 315.87 | 8.40 | 38.27 |
| 0.5 | 96 | 2948.20–3082.98 | 649.87–663.64 | 486.90 | 160.64 | 4.39 | 6.99 |
| 0.5 | 192 | 5766.79–6144.32 | 1390.44–1417.89 | 1070.09 | 326.14 | 8.73 | 41.67 |

Shared Metal library initialization was 37.03 ms, also outside warm solve time. CPU setup costs and all per-trial measurements are preserved in the CSV. This is a local run with ordinary system activity, not a thermally or scheduler-isolated measurement; CPU timings in particular show variation.

## RHS batching

For mass 0.1, the highest measured Metal throughput is 109.3 RHS/s at 192 RHS. Compare this with the timing spread above before interpreting a small difference between batch sizes.

For mass 0.5, the highest measured Metal throughput is 147.4 RHS/s at 96 RHS. Compare this with the timing spread above before interpreting a small difference between batch sizes.

These are independent batched solves, not block BiCGStab or recycling of a common Krylov basis. Additional RHS amortize control overhead and improve utilization, but increase vector traffic and FP64 CPU work. They do not guarantee higher throughput. The earlier SPD throughput numbers do not predict Wilson throughput or iteration counts.

## Validation and scope

All 36 measured complete solves met the per-RHS true-residual target. Every component was checked on the CPU in FP64 against the original FP64 Wilson operator. Known-solution relative error was at most 8.42e-10. The separate nine-system self-test checks gamma anticommutation, SU(3) unitarity and determinant, gamma5 Hermiticity, non-Hermiticity, direct Wilson versus encoded CPU/Metal products, free-field Fourier symbol, odd RHS, and a zero RHS. CPU and Metal reductions have different orders, so BiCGStab iteration counts can differ.

No clover term, even/odd Schur complement, multigrid, domain decomposition, physical negative bare mass, or production gauge ensemble is included. No Chroma/QDP binary was used as a numerical oracle. Operator normalization and the identity/gamma Kronecker decomposition follow the inspected local source, with an explicitly documented Euclidean spin basis. This is closer to mgproton.cc's workload, not a production substitute.

See [WILSON.md](WILSON.md) for equations, source references, implementation details, and reproduction.

Raw data: [results/wilson.csv](results/wilson.csv), [convergence traces](results/wilson.csv.trace.csv).
