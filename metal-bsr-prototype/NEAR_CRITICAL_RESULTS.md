# Clover-Wilson comparison near the estimated gap minimum

Tested bare mass **m₀=-0.12** on the same 16⁴ synthetic gauge field (strength 0.25), c_sw=1, QDP conventions, and antiperiodic time. The chosen mass is on the positive side of a finite-configuration low-gap minimum estimated near -0.14. Full residuals of the original reconstructed clover-Wilson equation are computed on the CPU in FP64 for every component and every RHS, with target 1e-10. Inner work and even/odd preparation/reconstruction are FP32 on the selected processor. Hardware: Apple M1 Max, 24 GPU cores, 32 GB, macOS 26.7.1. Built with Apple clang 21.0.0 using `-O3 -mcpu=apple-m1 -ffp-contract=fast`; full flags are in the Makefile.

## Complete solve times

Medians of 3 paired trials after a full resident-tile warmup per method. **96 RHS are resident; 192 RHS use two sequential 96-RHS tiles.** The first actual measured tile supplies the 96-RHS result and the sum of both distinct measured tiles supplies the 192-RHS result. Every tile starts from zero. Metal runs first, followed by CPU; method workspaces are allocated sequentially, retaining the same operator, hierarchy and common FP64 outer buffers.

| RHS | CPU BiCG s | CPU smoother FG s | CPU MG s | Metal BiCG s | Metal smoother FG s | Metal MG s |
|---:|---:|---:|---:|---:|---:|---:|
| 96 | 142.388 | 178.161 | 52.888 | 29.963 | 34.649 | 10.681 |
| 192 | 327.531 | 355.462 | 116.162 | 63.797 | 71.466 | 21.922 |

Solve time includes resets, all fine/coarse products, vector operations/reductions, smoothing, small CPU QR, GPU submission/waits, FP32 preparation and reconstruction, FP64 updates, and full CPU FP64 residual checks. Hierarchy/workspace/problem setup and manufactured RHS construction are separate. Initial common Metal library compilation is outside solve timing. Known-solution error checking follows residual acceptance outside timing. CPU uses eight workers; setup BLAS is single-threaded. These timings have ordinary local system activity and no thermal/scheduling isolation.

| RHS | CPU BiCG / MG | CPU smoother / MG | Metal BiCG / MG | Metal smoother / MG | MG CPU / Metal |
|---:|---:|---:|---:|---:|---:|
| 96 | 2.69× | 3.37× | 2.81× | 3.24× | 4.95× |
| 192 | 2.82× | 3.06× | 2.91× | 3.26× | 5.30× |

The first four ratios compare each baseline with MG; the last is the Metal speedup over CPU for MG. Medians and timing ranges describe these runs; small differences within the ranges need caution.

| RHS | Method | Min–max solve s | Full CPU FP64 residual s | FP32 inner s |
|---:|---|---:|---:|---:|
| 96 | cpu_bicgstab | 141.624–144.880 | 2.148 | 140.176 |
| 96 | cpu_smooth | 173.259–192.995 | 3.604 | 173.909 |
| 96 | cpu_sparse_mg | 52.570–68.193 | 2.397 | 50.295 |
| 96 | metal_bicgstab | 29.590–30.634 | 2.200 | 27.662 |
| 96 | metal_smooth | 34.286–37.778 | 2.321 | 32.272 |
| 96 | metal_sparse_mg | 10.461–10.928 | 2.416 | 8.028 |
| 192 | cpu_bicgstab | 322.365–328.370 | 4.260 | 323.077 |
| 192 | cpu_smooth | 339.341–365.940 | 6.620 | 347.251 |
| 192 | cpu_sparse_mg | 106.116–137.701 | 6.296 | 109.050 |
| 192 | metal_bicgstab | 63.307–65.032 | 4.384 | 59.181 |
| 192 | metal_smooth | 69.503–76.320 | 4.890 | 66.149 |
| 192 | metal_sparse_mg | 20.910–23.627 | 5.283 | 16.327 |

For Metal MG at 96 RHS, full CPU FP64 residual checks take 2.416 s, or 22.6% of complete solve time. This cost remains in every timing above.

## Convergence and actual work

The inner target is 0.002 for every method. Fine FGMRES restart remains 8, and smoother-only uses two clover-block Jacobi steps at weight 1.0. MG has one pre/post step at the same weight, 12 vectors per chirality, 4⁴ aggregates, and eight setup filtering sweeps at weight 0.65. Its sparse coarse solve uses a stronger 16-iteration cap and target preconditioned residual 0.1, chosen by the pilot below. BiCGStab retains its odd-Schur baseline, and Metal checks convergence every iteration (chunk=1). All methods have the same full FP64 target.

| RHS | Method | Outer corrections/tile | Summed Krylov steps/tile | Total fine S calls | Coarse SpMV | Worst full FP64 residual |
|---:|---|---:|---:|---:|---:|---:|
| 96 | cpu_bicgstab | 4 | 280 | 556 | 0 | 9.31e-12 |
| 96 | cpu_smooth | 4 | 249 | 531 | 0 | 9.39e-12 |
| 96 | cpu_sparse_mg | 4 | 30 | 94 | 552 | 1.01e-12 |
| 96 | metal_bicgstab | 4 | 270 | 540 | 0 | 9.28e-12 |
| 96 | metal_smooth | 4 | 249 | 531 | 0 | 9.39e-12 |
| 96 | metal_sparse_mg | 4 | 30 | 94 | 990 | 1.01e-12 |
| 192 | cpu_bicgstab | 4 | 280–358 | 1270 | 0 | 9.96e-12 |
| 192 | cpu_smooth | 4 | 247–249 | 1058 | 0 | 9.39e-12 |
| 192 | cpu_sparse_mg | 4 | 30 | 188 | 1121 | 1.01e-12 |
| 192 | metal_bicgstab | 4 | 270–313 | 1166 | 0 | 1e-11 |
| 192 | metal_smooth | 4 | 247–249 | 1058 | 0 | 9.39e-12 |
| 192 | metal_sparse_mg | 4 | 30 | 188 | 1980 | 1.01e-12 |

Krylov counts sum per-RHS maxima across outer corrections. BiCGStab and FGMRES iterations are different work. Every fine Schur application is counted, including smoothing, refreshes and GPU half-step dispatches. Calls at 192 RHS sum both 96-wide tiles. GPU coarse counts include masked/padded dispatches. Actual coarse residuals are diagnostic and do not accept a solve.

Example full FP64 residual history for global RHS offset 96, trial 0:

| Outer correction | CPU BiCG | CPU smoother FG | CPU MG | Metal BiCG | Metal smoother FG | Metal MG |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 0.00131 | 0.00131 | 0.000607 | 0.00131 | 0.00131 | 0.000607 |
| 2 | 2.55e-06 | 2.56e-06 | 6.25e-07 | 2.55e-06 | 2.56e-06 | 6.25e-07 |
| 3 | 5.08e-09 | 4.89e-09 | 7.57e-10 | 5.06e-09 | 4.89e-09 | 7.57e-10 |
| 4 | 9.96e-12 | 9.06e-12 | 9.87e-13 | 1e-11 | 9.06e-12 | 9.87e-13 |

### Comparison with the earlier m₀=+0.1 test

Historical positive-mass data provide a convergence/work reference. These are not additional fresh timing samples. The near-critical MG configuration also has a stronger coarse budget; that change is part of the reported configuration.

| Method | Old +0.1 fine S / 96 RHS | New -0.12 fine S / 96 RHS | Old +0.1 Krylov steps | New -0.12 Krylov steps |
|---|---:|---:|---:|---:|
| cpu_bicgstab | 143 | 556 | 73 | 280 |
| cpu_smooth | 132 | 531 | 61 | 249 |
| cpu_sparse_mg | 76 | 94 | 24 | 30 |
| metal_bicgstab | 146 | 540 | 73 | 270 |
| metal_smooth | 132 | 531 | 61 | 249 |
| metal_sparse_mg | 76 | 94 | 24 | 30 |

## Spectral location and its limits

The separate FP64 diagnostic applies `Q=gamma5*M` and estimates the low spectrum of `Q²=M†M` using fully reorthogonalized Lanczos. The gap convention is the minimum absolute Q eigenvalue ([Del Debbio et al., section 2](https://arxiv.org/abs/hep-lat/0512021)). Candidate vectors are reconstructed and their actual full-lattice FP64 Q² eigen-residual is evaluated. Free-field tests match the analytic Wilson dispersion.

| Bare mass | Steps | Start seed | Singular Rayleigh estimate | Actual Q² eigen-residual | Stationary-mass diagnostic | Reached 1e-5 eigen-residual? |
|---:|---:|---:|---:|---:|---:|---|
| 0.1 | 300 | 12371 | 0.2773299 | 0.00309 | -0.107411 | no |
| -0.1 | 300 | 12371 | 0.1820169 | 0.00147 | -0.120454 | no |
| -0.12 | 300 | 12371 | 0.1795421 | 0.000481 | -0.137860 | no |
| -0.14 | 300 | 12371 | 0.1785361 | 0.000455 | -0.139941 | no |
| -0.14 | 500 | 92871 | 0.1785218 | 0.000101 | -0.139826 | no |

The scan estimates a low-gap minimum near -0.14, with m=-0.12 on its positive side. Different starts/lengths near -0.14 give consistent approximate Rayleigh values, but the finite-volume runs do not reach their requested eigen-residual tolerance. Their values are Rayleigh upper estimates, so no precise critical mass or certified global spectral bound is claimed. The stationary-mass diagnostic is `m-Re<v,M(m)v>` for the normalized candidate; it locates the minimum of the mass quadratic for that candidate. This is a fixed-field spectral proxy, with antiperiodic time retained. A physical PCAC critical mass has not been determined. The finite synthetic field still has a substantial estimated gap near the minimum. **The full solve acceptance remains 1e-10 in FP64 regardless of these approximate spectral diagnostics.**

## Coarse solver screening

The four-step coarse cap from the positive-mass test was inadequate near this mass. The small-RHS screening uses the same field and strict full residual checks; it selects coarse accuracy and fine work, not a final timing result. Some screening diagnostics overlapped other functional/spectral checks. The final comparison runs methods sequentially.

| Coarse cap | Resident RHS | Fine Krylov steps | Fine S | Worst actual coarse residual | Full FP64 residual |
|---:|---:|---:|---:|---:|---:|
| 4 | 12 | 80 | 252 | 73.88 | 3.71e-12 |
| 8 | 12 | 30 | 94 | 1.399 | 1.29e-12 |
| 12 | 12 | 30 | 94 | 0.09968 | 9.92e-13 |
| 16 | 12 | 30 | 94 | 0.09968 | 9.92e-13 |

Caps 12 and 16 gave the same pilot convergence with a good actual coarse residual; 16 was chosen for the wider RHS batch. CPU stops submitting coarse products after convergence. Metal encodes the cap with per-RHS masking, so extra padded work is inside the measured MG time. FGMRES stores the preconditioned basis to permit the inexact/variable correction.

## Setup and many RHS

| RHS | MG hierarchy s | CPU hierarchy-only break-even vs best baseline | Metal hierarchy-only break-even vs best baseline |
|---:|---:|---:|---:|
| 96 | 2.524 | 0.03 batches | 0.13 batches |
| 192 | 2.524 | 0.01 batches | 0.06 batches |

These estimates divide hierarchy setup by observed warm batch time saved, excluding differences in method-specific workspace setup and first-use costs. The hierarchy is reused for all RHS, processors and trials at this mass; more RHS amortize it.

| Method | Inner workspace s | Common outer buffers s | Common problem setup s | Peak process RSS GiB | Worst coarse true residual |
|---|---:|---:|---:|---:|---:|
| cpu_bicgstab | 0.485 | 0.899 | 2.675 | 10.92 | 0 |
| cpu_smooth | 1.558 | 0.899 | 2.675 | 16.21 | 0 |
| cpu_sparse_mg | 0.685 | 0.899 | 2.675 | 16.68 | 0.101 |
| metal_bicgstab | 0.055 | 0.899 | 2.675 | 7.60 | 0 |
| metal_smooth | 0.102 | 0.899 | 2.675 | 7.30 | 0 |
| metal_sparse_mg | 0.117 | 0.899 | 2.675 | 7.60 | 0.101 |

RSS is the cumulative process high-water mark, not an isolated per-method measurement or total system memory. The sparse matrix retains 6,144 coarse unknowns and 37.125 MiB of FP32 coefficients; large fine vectors dominate memory.

## Validation and reproduction

All 36 paired logical records passed every-RHS full CPU FP64 residual acceptance. They derive from 36 physical 96-wide tile solves (3,456 RHS executions). Worst final relative residual: 1e-11; maximum known-solution error: 1.09e-10. CSV acceptance values were audited against final per-tile traces. Negative-mass small-lattice tests cover all six paths and zero RHS. The actual 16⁴ m=-0.12 sparse hierarchy passes full FP64 Galerkin/gamma5 checks (relative Galerkin error about 1.8e-15), plus strict CPU/Metal coarse solves. Local clover Cholesky factors are positive; the reported smallest even-factor pivot is about 3.4854.

```sh
make build/critical_probe build/metal_solver_compare
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \
  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass -.12 --samples 3 \
  --threads 8 --csw 1 --bi-chunk 1 --coarse-iters 16 --coarse-tol .1 \
  --metal-first --csv results/nearcritical-compare16.csv
python3 nearcritical_summary.py
```

See [NEAR_CRITICAL.md](NEAR_CRITICAL.md) for spectral methodology, pilot details and the timing protocol. These are synthetic-field results. Chroma and SuperbBLAS sources are unchanged.

Raw [measurements](results/nearcritical-compare16.csv), [full FP64 traces](results/nearcritical-compare16.csv.trace.csv), [run log](results/nearcritical-compare16.log), [spectral scan](results/critical-refine.csv), [independent spectral check](results/critical-independent.csv).
