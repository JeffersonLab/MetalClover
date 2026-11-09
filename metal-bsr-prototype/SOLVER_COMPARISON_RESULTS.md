# BiCGStab, smoother FGMRES, and sparse multigrid on 16⁴

Fresh comparison on the Apple M1 Max (24 GPU cores, 32 GB, macOS 26.7.1). All methods use the same 16⁴ synthetic clover-Wilson operator, c_sw=1, gauge strength 0.25, QDP spin conventions, antiperiodic time, and positive bare masses 0.1/0.5. All RHS, odd Schur reduction/reconstruction, and stopping targets agree across methods. Inner work is FP32; accumulation and every component of the full reconstructed residual are CPU FP64. Every RHS must meet a true relative residual of 1e-10.

## Complete solve times

Medians of three trials from zero, after one complete resident-tile warmup per method. **96 RHS are resident. 192 RHS use two sequential 96-RHS tiles.** Within a trial the first measured tile supplies the 96-RHS time, and both actual measured tiles supply the 192-RHS total. Thus the two sizes are paired prefix measurements; 192 is not an extrapolation and is not 192 simultaneously resident RHS.

| Mass | RHS | CPU BiCG s | CPU smoother FG s | CPU sparse MG s | Metal BiCG s | Metal smoother FG s | Metal sparse MG s |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 46.517 | 57.962 | 44.275 | 9.990 | 11.039 | 7.821 |
| 0.1 | 192 | 92.975 | 109.663 | 81.920 | 19.724 | 21.971 | 15.722 |
| 0.5 | 96 | 19.806 | 32.916 | 23.467 | 5.941 | 7.180 | 7.808 |
| 0.5 | 192 | 39.369 | 57.590 | 49.203 | 11.624 | 14.260 | 17.726 |

Timings include zero/reset work, FP32 preparation/reconstruction, all fine and coarse products, vector operations/reductions, small CPU QR, GPU submissions/waits, FP64 solution updates and full FP64 residual checks. Setup and manufactured RHS construction are outside solve time and reported separately. Known-solution error is checked after acceptance outside timing. CPU uses eight workers. Methods run sequentially, rotating order across mass/processor cases, with processor order reversed between masses. Ordinary local system activity was present; thermal and scheduling state were not isolated.

| Mass | RHS | Fastest CPU median | Fastest Metal median | CPU BiCG / MG | CPU smoother / MG | Metal BiCG / MG | Metal smoother / MG |
|---:|---:|---|---|---:|---:|---:|---:|
| 0.1 | 96 | sparse MG | sparse MG | 1.05× | 1.31× | 1.28× | 1.41× |
| 0.1 | 192 | sparse MG | sparse MG | 1.13× | 1.34× | 1.25× | 1.40× |
| 0.5 | 96 | BiCGStab | BiCGStab | 0.84× | 1.40× | 0.76× | 0.92× |
| 0.5 | 192 | BiCGStab | BiCGStab | 0.80× | 1.17× | 0.66× | 0.80× |

A baseline/MG ratio above one favors sparse MG; below one favors the baseline. Fastest refers to the observed median. Small differences within the timing spread do not establish a robust ranking.

| Mass | RHS | BiCG CPU / Metal | Smoother FG CPU / Metal | Sparse MG CPU / Metal |
|---:|---:|---:|---:|---:|
| 0.1 | 96 | 4.66× | 5.25× | 5.66× |
| 0.1 | 192 | 4.71× | 4.99× | 5.21× |
| 0.5 | 96 | 3.33× | 4.58× | 3.01× |
| 0.5 | 192 | 3.39× | 4.04× | 2.78× |

## Convergence and actual operator work

BiCGStab and FGMRES iterations represent different work. Fine S counts include smoothing, Arnoldi products, true FP32 refreshes, and every dispatched BiCGStab product. Each call applies to one resident 96-RHS tile; 192-RHS rows sum both tiles. GPU BiCGStab checks convergence every iteration (chunk=1), minimizing iteration padding; its second half-step is still dispatched when RHS are already converged. Krylov counts sum the maximum per-RHS iteration count over outer corrections; ranges cover resident tiles/trials. All methods retain full FP64 outer checks.

| Mass | RHS | Method | Outer corrections/tile | Krylov steps/tile | Total fine S | Coarse SpMV | Worst full FP64 residual |
|---:|---:|---|---:|---:|---:|---:|---:|
| 0.1 | 96 | cpu_bicgstab | 4 | 73 | 143 | 0 | 8.26e-12 |
| 0.1 | 96 | cpu_smooth | 4 | 61 | 132 | 0 | 6.51e-12 |
| 0.1 | 96 | cpu_sparse_mg | 4 | 24 | 76 | 135 | 9.59e-13 |
| 0.1 | 96 | metal_bicgstab | 4 | 73 | 146 | 0 | 7.88e-12 |
| 0.1 | 96 | metal_smooth | 4 | 61 | 132 | 0 | 6.51e-12 |
| 0.1 | 96 | metal_sparse_mg | 4 | 24 | 76 | 216 | 9.59e-13 |
| 0.1 | 192 | cpu_bicgstab | 4 | 70–73 | 279 | 0 | 8.26e-12 |
| 0.1 | 192 | cpu_smooth | 4 | 61 | 264 | 0 | 6.51e-12 |
| 0.1 | 192 | cpu_sparse_mg | 4 | 24 | 152 | 269 | 9.68e-13 |
| 0.1 | 192 | metal_bicgstab | 4 | 70–73 | 286 | 0 | 8.41e-12 |
| 0.1 | 192 | metal_smooth | 4 | 61 | 264 | 0 | 6.51e-12 |
| 0.1 | 192 | metal_sparse_mg | 4 | 24 | 152 | 432 | 9.68e-13 |
| 0.5 | 96 | cpu_bicgstab | 4 | 31 | 59 | 0 | 7.6e-12 |
| 0.5 | 96 | cpu_smooth | 4 | 27 | 58 | 0 | 2.1e-12 |
| 0.5 | 96 | cpu_sparse_mg | 4 | 19 | 61 | 59 | 2.38e-13 |
| 0.5 | 96 | metal_bicgstab | 4 | 31 | 62 | 0 | 7.6e-12 |
| 0.5 | 96 | metal_smooth | 4 | 27 | 58 | 0 | 2.1e-12 |
| 0.5 | 96 | metal_sparse_mg | 4 | 19 | 61 | 171 | 2.38e-13 |
| 0.5 | 192 | cpu_bicgstab | 4 | 29–31 | 115 | 0 | 7.74e-12 |
| 0.5 | 192 | cpu_smooth | 4 | 27 | 116 | 0 | 2.1e-12 |
| 0.5 | 192 | cpu_sparse_mg | 4 | 19 | 122 | 118 | 2.39e-13 |
| 0.5 | 192 | metal_bicgstab | 4 | 29–31 | 120 | 0 | 7.73e-12 |
| 0.5 | 192 | metal_smooth | 4 | 27 | 116 | 0 | 2.1e-12 |
| 0.5 | 192 | metal_sparse_mg | 4 | 19 | 122 | 342 | 2.39e-13 |

Smoother FGMRES uses two clover-block Jacobi steps at weight 1.0, with restart 8. Sparse MG uses one pre/post step at the same weight and restart, 12 basis vectors per chirality, 4⁴ aggregates, and eight FP64 filtering sweeps at setup weight 0.65. Its 6,144-unknown coarse operator has 8,448 dense 24×24 sparse blocks. Coarse BiCGStab has local block preconditioning, target 0.1 in the preconditioned residual and a four-iteration cap. Flexible GMRES permits this inexact, variable correction. The actual unpreconditioned FP32 coarse residual is recomputed and reported in raw data; it may exceed 0.1. It never replaces the full FP64 stopping test. GPU coarse counts include masked/padded dispatches.

Full CPU FP64 history for the second 96-RHS tile (global offset 96), mass 0.1, trial 0:

| Outer correction | CPU BiCG | CPU smoother FG | CPU MG | Metal BiCG | Metal smoother FG | Metal MG |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.00132 | 0.000999 | 0.000514 | 0.00132 | 0.000999 | 0.000514 |
| 2 | 2.57e-06 | 1.96e-06 | 8.65e-07 | 2.58e-06 | 1.96e-06 | 8.65e-07 |
| 3 | 4.71e-09 | 3.34e-09 | 1.68e-09 | 4.82e-09 | 3.34e-09 | 1.68e-09 |
| 4 | 8.21e-12 | 6.46e-12 | 9.68e-13 | 8.41e-12 | 6.46e-12 | 9.68e-13 |

## Setup and an unlimited supply of RHS

Each mass-specific hierarchy is built once and reused across all methods, processors and trials. Baselines need no hierarchy; their hierarchy setup field is zero even though the shared experiment also prepares MG. Fine-operator/problem generation and FP64 outer storage are common. Each solver allocates its own inner workspace sequentially. Allocation timing does not include all first-use page/pipeline costs, which warmup covers.

| Mass | RHS | MG hierarchy s | CPU best baseline | CPU hierarchy-only break-even batches | Metal best baseline | Metal hierarchy-only break-even batches |
|---:|---:|---:|---|---:|---|---:|
| 0.1 | 96 | 3.028 | BiCGStab | 1.4 | BiCGStab | 1.4 |
| 0.1 | 192 | 3.028 | BiCGStab | 0.3 | BiCGStab | 0.8 |
| 0.5 | 96 | 5.722 | BiCGStab | No warm gain | BiCGStab | No warm gain |
| 0.5 | 192 | 5.722 | BiCGStab | No warm gain | BiCGStab | No warm gain |

The estimate divides hierarchy setup by measured warm batch time saved relative to the faster baseline, excluding method-specific workspace/compilation differences. It is sensitive to small timing differences. More RHS can amortize setup when MG has a warm-time gain; it cannot remove a per-batch slowdown.

| Mass | Method | Inner workspace setup s | Common outer setup s | Common problem setup s | Peak process RSS GiB |
|---:|---|---:|---:|---:|---:|
| 0.1 | cpu_bicgstab | 0.620 | 0.942 | 3.455 | 10.90 |
| 0.1 | cpu_smooth | 2.174 | 0.942 | 3.455 | 16.21 |
| 0.1 | cpu_sparse_mg | 0.585 | 0.942 | 3.455 | 16.68 |
| 0.1 | metal_bicgstab | 0.059 | 0.942 | 3.455 | 16.68 |
| 0.1 | metal_smooth | 0.566 | 0.942 | 3.455 | 16.68 |
| 0.1 | metal_sparse_mg | 0.165 | 0.942 | 3.455 | 16.68 |
| 0.5 | cpu_bicgstab | 0.066 | 0.037 | 2.381 | 16.68 |
| 0.5 | cpu_smooth | 0.839 | 0.037 | 2.381 | 16.68 |
| 0.5 | cpu_sparse_mg | 3.690 | 0.037 | 2.381 | 16.68 |
| 0.5 | metal_bicgstab | 0.072 | 0.037 | 2.381 | 16.68 |
| 0.5 | metal_smooth | 0.086 | 0.037 | 2.381 | 16.68 |
| 0.5 | metal_sparse_mg | 0.060 | 0.037 | 2.381 | 16.68 |

Process RSS is the cumulative high-water mark within the run, rather than an isolated per-method memory measurement. The sparse FP32 coefficient matrix occupies 37.125 MiB, with a separate FP64 copy and basis/local inverses. Large fine Krylov vectors dominate memory, hence sequential method allocation and 96-wide tiling.

## Timing spread and FP64 cost

| Mass | RHS | Method | Min–max solve s | Median full FP64 residual s | Median FP32 inner s |
|---:|---:|---|---:|---:|---:|
| 0.1 | 96 | cpu_bicgstab | 45.840–48.484 | 3.425 | 43.303 |
| 0.1 | 96 | cpu_smooth | 53.948–60.589 | 3.906 | 53.380 |
| 0.1 | 96 | cpu_sparse_mg | 33.001–44.323 | 3.809 | 39.765 |
| 0.1 | 96 | metal_bicgstab | 9.954–9.999 | 2.224 | 7.642 |
| 0.1 | 96 | metal_smooth | 10.795–11.268 | 2.436 | 8.468 |
| 0.1 | 96 | metal_sparse_mg | 7.775–7.961 | 2.404 | 5.275 |
| 0.1 | 192 | cpu_bicgstab | 91.046–96.203 | 6.700 | 85.457 |
| 0.1 | 192 | cpu_smooth | 106.525–117.475 | 7.923 | 99.230 |
| 0.1 | 192 | cpu_sparse_mg | 65.540–94.167 | 7.485 | 73.349 |
| 0.1 | 192 | metal_bicgstab | 19.633–19.762 | 4.482 | 14.992 |
| 0.1 | 192 | metal_smooth | 21.784–22.273 | 4.769 | 16.934 |
| 0.1 | 192 | metal_sparse_mg | 15.523–16.119 | 4.804 | 10.648 |
| 0.5 | 96 | cpu_bicgstab | 19.656–20.448 | 2.312 | 17.352 |
| 0.5 | 96 | cpu_smooth | 26.262–34.122 | 4.115 | 27.921 |
| 0.5 | 96 | cpu_sparse_mg | 23.163–24.582 | 2.258 | 21.079 |
| 0.5 | 96 | metal_bicgstab | 5.828–5.949 | 2.310 | 3.490 |
| 0.5 | 96 | metal_smooth | 6.538–7.587 | 2.365 | 4.620 |
| 0.5 | 96 | metal_sparse_mg | 7.294–10.614 | 2.489 | 5.222 |
| 0.5 | 192 | cpu_bicgstab | 39.212–39.537 | 4.641 | 34.443 |
| 0.5 | 192 | cpu_smooth | 51.578–63.309 | 6.938 | 49.460 |
| 0.5 | 192 | cpu_sparse_mg | 46.040–56.064 | 4.547 | 44.427 |
| 0.5 | 192 | metal_bicgstab | 11.562–11.892 | 4.577 | 6.804 |
| 0.5 | 192 | metal_smooth | 13.986–21.883 | 5.352 | 8.725 |
| 0.5 | 192 | metal_sparse_mg | 15.726–19.551 | 5.476 | 11.442 |

## Validation and reproduction

All 72 paired logical batch records passed every-RHS full CPU FP64 acceptance. They derive from 72 physically measured 96-RHS tile solves (6,912 RHS executions), with shared first tiles between the reported logical sizes. Worst final relative residual: 8.41e-12; maximum known-solution error: 5.86e-11. Every CSV acceptance value was audited against the final per-tile FP64 trace values. The 48-record small comparison self-test passed all six paths, two volumes and masses, and zero RHS. The 40-solve established clover regression passed after ensuring GPU BiCGStab allocates a full direction workspace when optional multiply inputs have been released. Sparse Galerkin validation from the previous study remains applicable: no hierarchy or fine/coarse operator algorithm changed.

```sh
make build/metal_solver_compare
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare --self-test
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \
  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass .1,.5 --csw 1 \
  --samples 3 --threads 8 --bi-chunk 1 --csv results/solver-compare16.csv
python3 solver_compare_summary.py
```

See [SOLVER_COMPARISON.md](SOLVER_COMPARISON.md) for algorithms and the paired timing protocol. Results describe synthetic positive-mass fields on this M1 Max, not a production ensemble or critical-mass scaling study. Chroma and SuperbBLAS sources are unchanged.

Raw [measurements](results/solver-compare16.csv), [full FP64 traces](results/solver-compare16.csv.trace.csv), [run log](results/solver-compare16.log).
