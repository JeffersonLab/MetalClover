# Even/odd Wilson benchmark on M1 Max

Measured locally on November 3, 2026: M1 Max, 24 GPU cores, 32 GB, macOS 26.7.1, AC power, eight CPU worker tasks. Each entry is the median of three complete solves from zero, after one full warmup per method. All six methods are measured in the same run, with rotating order. The links and RHS are identical across methods within each case.

The retained operator is the bare Wilson matrix with site-dependent synthetic SU(3) links (strength 0.25), periodic space, and antiperiodic time. Even/odd uses packed checkerboards and the normalized Schur complement `S = I - H_eo H_oe / (m+4)^2`. Each Schur product uses two half-volume hopping passes. All paths use complex BiCGStab directly, with no normal equations.

**Every RHS is accepted only after full odd-site reconstruction and a CPU FP64 check of `||b-D64*x||/||b|| <= 1e-10`.** Mixed paths use FP32 inner tolerance 0.002 and FP64 correction accumulation. EO RHS preparation and odd reconstruction are also CPU FP64 work, included in complete solve time. Initial allocation, operator upload, and compilation are outside warm solve time.

## Complete runtime: with and without even/odd

CPU best means the faster measured CPU FP64 or CPU mixed median for that preconditioner. EO gain means unpreconditioned time divided by preconditioned time; a value below one is a slowdown.

| Mass | Lattice | RHS | CPU best, none (ms) | CPU best, EO (ms) | Metal, none (ms) | Metal, EO (ms) | Metal EO gain | CPU EO / Metal EO | EO Metal ms/RHS |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 8⁴ | 96 | 4712.66 | 2607.91 | 1006.71 | 788.28 | 1.28× | 3.31× | 8.211 |
| 0.1 | 8⁴ | 192 | 9059.75 | 4903.92 | 2046.48 | 1734.81 | 1.18× | 2.83× | 9.035 |
| 0.5 | 8⁴ | 96 | 3235.01 | 1749.91 | 709.02 | 752.83 | 0.94× | 2.32× | 7.842 |
| 0.5 | 8⁴ | 192 | 5949.62 | 3439.42 | 1483.99 | 1355.46 | 1.09× | 2.54× | 7.060 |

## All CPU baselines and convergence

Counts are `outer corrections / summed maximum per-RHS inner iterations / dispatched iterations`. A Schur iteration is not a half-cost Wilson iteration: each operator application makes two half-volume hopping passes. BiCGStab needs up to two operator applications per iteration. Metal pads to four-iteration command chunks; the CPU can stop early. Timing includes that work.

| Mass | RHS | Preconditioner | Backend | Total (ms) | Corrections / inner / dispatched | Worst full residual |
|---:|---:|---|---|---:|---:|---:|
| 0.1 | 96 | none | cpu_fp64 | 7814.46 | 1 / 72 / 72 | 9.99e-11 |
| 0.1 | 96 | none | cpu_mixed | 4712.66 | 4 / 86 / 86 | 1.3e-11 |
| 0.1 | 96 | none | metal_mixed | 1006.71 | 4 / 85 / 92 | 1.31e-11 |
| 0.1 | 96 | even_odd | cpu_fp64 | 3580.42 | 1 / 36 / 36 | 6.58e-11 |
| 0.1 | 96 | even_odd | cpu_mixed | 2607.91 | 4 / 44 / 44 | 7.64e-12 |
| 0.1 | 96 | even_odd | metal_mixed | 788.28 | 4 / 44 / 48 | 7.64e-12 |
| 0.1 | 192 | none | cpu_fp64 | 15084.43 | 1 / 71 / 71 | 9.99e-11 |
| 0.1 | 192 | none | cpu_mixed | 9059.75 | 4 / 87 / 87 | 1.44e-11 |
| 0.1 | 192 | none | metal_mixed | 2046.48 | 4 / 87 / 92 | 1.34e-11 |
| 0.1 | 192 | even_odd | cpu_fp64 | 7166.07 | 1 / 36 / 36 | 6.58e-11 |
| 0.1 | 192 | even_odd | cpu_mixed | 4903.92 | 4 / 46 / 46 | 8.42e-12 |
| 0.1 | 192 | even_odd | metal_mixed | 1734.81 | 4 / 46 / 52 | 8.4e-12 |
| 0.5 | 96 | none | cpu_fp64 | 4090.04 | 1 / 37 / 37 | 9.8e-11 |
| 0.5 | 96 | none | cpu_mixed | 3235.01 | 4 / 59 / 59 | 1.25e-11 |
| 0.5 | 96 | none | metal_mixed | 709.02 | 4 / 58 / 64 | 1.26e-11 |
| 0.5 | 96 | even_odd | cpu_fp64 | 1749.91 | 1 / 19 / 19 | 6.2e-11 |
| 0.5 | 96 | even_odd | cpu_mixed | 1835.12 | 4 / 32 / 32 | 7.65e-12 |
| 0.5 | 96 | even_odd | metal_mixed | 752.83 | 4 / 32 / 44 | 7.65e-12 |
| 0.5 | 192 | none | cpu_fp64 | 7110.61 | 1 / 38 / 38 | 9.5e-11 |
| 0.5 | 192 | none | cpu_mixed | 5949.62 | 4 / 61 / 61 | 1.23e-11 |
| 0.5 | 192 | none | metal_mixed | 1483.99 | 4 / 61 / 72 | 1.32e-11 |
| 0.5 | 192 | even_odd | cpu_fp64 | 3439.42 | 1 / 19 / 19 | 6.45e-11 |
| 0.5 | 192 | even_odd | cpu_mixed | 3631.80 | 4 / 33 / 33 | 9.82e-11 |
| 0.5 | 192 | even_odd | metal_mixed | 1355.46 | 4 / 33 / 44 | 9.83e-11 |

## Where the preconditioned Metal time goes

The inner column excludes RHS preparation and reconstruction, which the raw `inner_ms` field includes. Component medians need not sum to the median total. Setup is the original allocation and checkerboard extraction/upload cost, not a cold-start solve.

| Mass | RHS | FP32 Schur BiCGStab (ms) | FP64 RHS reduction (ms) | FP64 odd reconstruction (ms) | Full FP64 residual checks (ms) | FP64 update (ms) | EO setup (ms) |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 405.15 | 96.93 | 95.68 | 182.23 | 6.22 | 67.87 |
| 0.1 | 192 | 850.64 | 242.17 | 185.79 | 365.04 | 11.74 | 143.01 |
| 0.5 | 96 | 368.90 | 96.13 | 97.95 | 182.29 | 6.33 | 99.74 |
| 0.5 | 192 | 609.18 | 180.02 | 171.83 | 320.04 | 11.24 | 90.30 |

## Full-Wilson convergence example: mass 0.1, 192 RHS

These are full-system residuals after reconstruction, not reduced-system residuals.

| Correction | Metal none | Metal even/odd |
|---:|---:|---:|
| 1 | 0.002 | 0.00131 |
| 2 | 3.98e-06 | 2.44e-06 |
| 3 | 7.46e-09 | 4.53e-09 |
| 4 | 1.34e-11 | 8.4e-12 |

## Timing spread and batching

| Mass | RHS | EO CPU mixed min–max (ms) | EO Metal min–max (ms) | EO Metal RHS/s |
|---:|---:|---:|---:|---:|
| 0.1 | 96 | 2421.81–2611.04 | 773.28–1097.83 | 121.8 |
| 0.1 | 192 | 4837.36–5107.00 | 1573.62–1769.57 | 110.7 |
| 0.5 | 96 | 1779.91–2449.59 | 749.10–889.23 | 127.5 |
| 0.5 | 192 | 3359.66–3641.47 | 1250.07–1372.24 | 141.6 |

Compare milliseconds per RHS and the timing spread when choosing 96 versus 192 RHS. All solves are independent batched BiCGStab; there is no shared Krylov basis or block solver. This was an ordinary local run, not a thermally or scheduler-isolated measurement. The CSV retains all individual trial times.

## Validation and scope

All 72 measured solves passed full FP64 residual and known-solution checks. Maximum known-solution relative error was 8.42e-10. The separate six-system self-test checks direct full-Wilson composition versus CPU/Metal Schur products, transformed RHS, odd reconstruction, odd RHS, zero RHS locking, duplicate neighbors, periodic and antiperiodic boundaries, and rejection of odd periodic extents. A separate free-field periodic case passed all six solver paths, and the nine-system unpreconditioned Wilson regression test passed after the shared-code changes.

The operator still has synthetic SU(3) links, positive bare mass, and no clover, multigrid, domain decomposition, or production gauge ensemble. Generic color/spin kernels are used on both processors; neither is a tuned half-spinor Dslash implementation. FP64 RHS/reconstruction on the CPU is a deliberate implementation choice and can limit net Metal speedup despite fewer inner iterations.

See [EVENODD.md](EVENODD.md) for formulas and reproduction.

Raw measurements: [results/evenodd.csv](results/evenodd.csv) and [convergence traces](results/evenodd.csv.trace.csv).
