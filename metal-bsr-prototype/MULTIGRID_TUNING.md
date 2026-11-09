# Clover multigrid tuning on M1 Max

Measured locally on November 4, 2026, on an Apple M1 Max, 24 GPU cores, 32 GB, macOS 26.7.1. This extends only the standalone prototype. The fine operator remains the Chroma-convention odd clover Schur complement with c_sw=1, positive bare masses 0.1/0.5, synthetic gauge strength 0.25, an 8⁴ lattice, and antiperiodic time. All CPU/Metal inner arithmetic and elimination are FP32. Every full residual component, including reconstructed even sites and clover terms, is evaluated on the CPU in FP64; every RHS must meet a relative residual of 1e-10.

## Selected settings

The selected hierarchy doubles the local basis from 6 to **12 vectors per chirality**, retains **4×4×4×4 aggregates**, and changes the cycle's clover-block Jacobi weight from **0.65 to 1.0**. One pre- and one post-smoothing step remain. The coarse matrix grows from 192² to 384². Eight FP64 test-vector filtering sweeps at the independent setup weight 0.65 suffice: longer filtering did not reduce the measured large-batch iteration counts. FGMRES restart stays 8, the inner target stays 0.002, and full FP64 acceptance stays 1e-10.

## Fine work and convergence

Counts are summed maximum per-RHS Krylov iterations across outer corrections, and actual batched fine Schur applications. CPU and Metal use identical counts for these multigrid cases. One Schur application includes both half-volume hopping operations and the local clover blocks. Counts include smoothing, Arnoldi work, and residual refreshes.

| Mass | RHS | Old → tuned Krylov iterations | Old → tuned fine S calls | Fine work reduction | Outer FP64 corrections |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 45 → 23 | 142 → 73 | 48.6% | 4 |
| 0.1 | 192 | 45 → 23 | 142 → 73 | 48.6% | 4 |
| 0.5 | 96 | 29 → 19 | 91 → 61 | 33.0% | 4 |
| 0.5 | 192 | 29 → 19 | 91 → 61 | 33.0% | 4 |

## Fresh complete solve timings

Medians of three measured solves from zero, following a full warmup for every method. Method order rotates within trials. Original and tuned configuration order alternates across mass/RHS cases; they run sequentially, never as competing benchmarks. CPU uses eight GCD workers and single-threaded Accelerate for coarse GEMM. Times include FP32 preparation/reconstruction, vector/reduction work, CPU small QR, GPU submissions/waits, FP64 accumulation, and full FP64 residual checks. Hierarchy/workspace/problem setup is excluded and reported separately. Ordinary local system activity was present.

| Mass | RHS | CPU old MG ms | CPU tuned MG ms | Metal old MG ms | Metal tuned MG ms | CPU tuning gain | Metal tuning gain | Tuned CPU / Metal |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 5661.96 | 2214.19 | 771.38 | 520.61 | 2.56× | 1.48× | 4.25× |
| 0.1 | 192 | 7786.88 | 4220.38 | 1534.16 | 1020.30 | 1.85× | 1.50× | 4.14× |
| 0.5 | 96 | 2550.37 | 1764.29 | 614.03 | 469.06 | 1.45× | 1.31× | 3.76× |
| 0.5 | 192 | 5230.25 | 3642.78 | 1175.09 | 985.64 | 1.44× | 1.19× | 3.70× |

### BiCGStab and the smoother-only control

BiCGStab is freshly measured with the original FP32 CPU/device elimination paths. It does not use the hierarchy. The smoother-only FGMRES control uses exactly the tuned two-step smoother at weight 1.0, with no coarse correction. This distinguishes stronger smoothing from the coarse-grid contribution.

| Mass | RHS | CPU BiCG ms | CPU smooth FG ms | CPU tuned MG ms | Metal BiCG ms | Metal smooth FG ms | Metal tuned MG ms |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 2454.65 | 2743.35 | 2214.19 | 604.72 | 615.32 | 520.61 |
| 0.1 | 192 | 4719.71 | 5402.85 | 4220.38 | 1225.99 | 1263.99 | 1020.30 |
| 0.5 | 96 | 1742.76 | 1675.82 | 1764.29 | 504.42 | 452.99 | 469.06 |
| 0.5 | 192 | 3439.91 | 3480.74 | 3642.78 | 999.64 | 911.89 | 985.64 |

A gain below 1 means the tuned two-level method is slower than that baseline. More available RHS can amortize setup only when a warm-time gain actually exists.

| Mass | RHS | CPU BiCG / tuned MG | Metal BiCG / tuned MG | CPU smooth FG / tuned MG | Metal smooth FG / tuned MG |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 1.11× | 1.16× | 1.24× | 1.18× |
| 0.1 | 192 | 1.12× | 1.20× | 1.28× | 1.24× |
| 0.5 | 96 | 0.99× | 1.08× | 0.95× | 0.97× |
| 0.5 | 192 | 0.94× | 1.01× | 0.96× | 0.93× |

## Full counts and numerical checks

Counts below are `outer corrections / Krylov iterations / fine S applications / coarse corrections`. BiCGStab and FGMRES iteration counts represent different work; GPU BiCGStab includes dispatch padding. FGMRES refreshes its actual FP32 Schur residual on restart, and the outer solver always checks the full FP64 equation.

| Mass | RHS | Setting | Method | Counts | Worst full FP64 residual |
|---:|---:|---|---|---:|---:|
| 0.1 | 96 | original | cpu_bicgstab | 4/46/90/0 | 8.27e-12 |
| 0.1 | 96 | original | metal_bicgstab | 4/46/104/0 | 8.27e-12 |
| 0.1 | 96 | original | cpu_fgmres_mg | 4/45/142/45 | 7.14e-12 |
| 0.1 | 96 | original | metal_fgmres_mg | 4/45/142/45 | 7.14e-12 |
| 0.1 | 96 | tuned | cpu_fgmres_smooth | 4/42/91/0 | 5.75e-12 |
| 0.1 | 96 | tuned | metal_fgmres_smooth | 4/42/91/0 | 5.75e-12 |
| 0.1 | 96 | tuned | cpu_fgmres_mg | 4/23/73/23 | 1.41e-12 |
| 0.1 | 96 | tuned | metal_fgmres_mg | 4/23/73/23 | 1.41e-12 |
| 0.1 | 192 | original | cpu_bicgstab | 4/47/90/0 | 8.7e-12 |
| 0.1 | 192 | original | metal_bicgstab | 4/47/112/0 | 8.71e-12 |
| 0.1 | 192 | original | cpu_fgmres_mg | 4/45/142/45 | 6.64e-12 |
| 0.1 | 192 | original | metal_fgmres_mg | 4/45/142/45 | 6.64e-12 |
| 0.1 | 192 | tuned | cpu_fgmres_smooth | 4/43/93/0 | 5.77e-12 |
| 0.1 | 192 | tuned | metal_fgmres_smooth | 4/43/93/0 | 5.77e-12 |
| 0.1 | 192 | tuned | cpu_fgmres_mg | 4/23/73/23 | 1.42e-12 |
| 0.1 | 192 | tuned | metal_fgmres_mg | 4/23/73/23 | 1.42e-12 |
| 0.5 | 96 | original | cpu_bicgstab | 4/33/63/0 | 5.31e-11 |
| 0.5 | 96 | original | metal_bicgstab | 4/33/88/0 | 5.31e-11 |
| 0.5 | 96 | original | cpu_fgmres_mg | 4/29/91/29 | 1.8e-12 |
| 0.5 | 96 | original | metal_fgmres_mg | 4/29/91/29 | 1.8e-12 |
| 0.5 | 96 | tuned | cpu_fgmres_smooth | 4/25/54/0 | 1.89e-12 |
| 0.5 | 96 | tuned | metal_fgmres_smooth | 4/25/54/0 | 1.89e-12 |
| 0.5 | 96 | tuned | cpu_fgmres_mg | 4/19/61/19 | 1.12e-13 |
| 0.5 | 96 | tuned | metal_fgmres_mg | 4/19/61/19 | 1.12e-13 |
| 0.5 | 192 | original | cpu_bicgstab | 4/33/63/0 | 6.7e-12 |
| 0.5 | 192 | original | metal_bicgstab | 4/33/88/0 | 6.7e-12 |
| 0.5 | 192 | original | cpu_fgmres_mg | 4/29/91/29 | 1.68e-12 |
| 0.5 | 192 | original | metal_fgmres_mg | 4/29/91/29 | 1.68e-12 |
| 0.5 | 192 | tuned | cpu_fgmres_smooth | 4/25/54/0 | 1.87e-12 |
| 0.5 | 192 | tuned | metal_fgmres_smooth | 4/25/54/0 | 1.87e-12 |
| 0.5 | 192 | tuned | cpu_fgmres_mg | 4/19/61/19 | 1.12e-13 |
| 0.5 | 192 | tuned | metal_fgmres_mg | 4/19/61/19 | 1.12e-13 |

Example full-system FP64 residual history, mass 0.1, 192 RHS, Metal, trial 0:

| Correction | Original MG | Tuned MG | Tuned smoother FG | BiCGStab |
|---:|---:|---:|---:|---:|
| 1 | 0.0012 | 0.00132 | 0.00119 | 0.0013 |
| 2 | 2.04e-06 | 1.24e-06 | 2e-06 | 2.56e-06 |
| 3 | 3.84e-09 | 1.29e-09 | 3.42e-09 | 4.53e-09 |
| 4 | 6.64e-12 | 1.42e-12 | 5.77e-12 | 8.71e-12 |

## Setup, throughput, and amortization

Hierarchy setup is CPU FP64 and independent of the benchmark RHS and known solution. Each final case rebuilds it to alternate the configuration order, then reuses it for all warmups/trials/methods. The large-batch finalist runs additionally reuse one mass-specific hierarchy across 96 and 192 RHS. It can be retained for subsequent batches on the same gauge, mass, clover coefficient, and geometry.

| Mass | RHS | Old hierarchy ms | Tuned hierarchy ms | Tuned Metal ms/RHS | Tuned Metal RHS/s | Tuned Metal FP64 residual ms | Hierarchy-only break-even batches vs Metal BiCG |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 618.54 | 1205.92 | 5.423 | 184.4 | 182.50 | 14.3 |
| 0.1 | 192 | 580.80 | 1319.73 | 5.314 | 188.2 | 350.77 | 6.4 |
| 0.5 | 96 | 603.26 | 1325.54 | 4.886 | 204.7 | 180.19 | 37.5 |
| 0.5 | 192 | 530.17 | 1204.16 | 5.134 | 194.8 | 397.82 | 86.0 |

Break-even is hierarchy setup divided by measured warm time saved per batch, without compilation/workspace setup; it is an estimate, particularly sensitive when the warm-time difference is small. CPU FP64 residual time remains inside all warm totals. Problem/workspace setup and per-trial timing components are in the raw CSV.

| Mass | RHS | CPU old MG min–max ms | CPU tuned MG min–max ms | Metal old MG min–max ms | Metal tuned MG min–max ms |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 4004.76–6361.36 | 2171.05–2268.03 | 758.71–899.25 | 504.98–531.60 |
| 0.1 | 192 | 7480.52–7825.86 | 4206.43–4272.18 | 1512.47–2172.69 | 1017.31–1045.57 |
| 0.5 | 96 | 2513.15–2632.38 | 1751.43–1830.21 | 565.65–676.32 | 442.66–471.76 |
| 0.5 | 192 | 4956.69–5547.65 | 3616.25–4115.11 | 1137.48–1309.87 | 938.31–1101.85 |

The old CPU mass-0.1/96-RHS timing has substantial spread; its 2.56× median tuning gain should be read with that variability. Fine Schur counts and numerical convergence are reproducible. At mass 0.5 the tuned smoother-only method remains faster than MG. The small Metal MG advantage over BiCGStab at 192 RHS/mass 0.5 is within the observed timing spread and does not establish a reliable setup break-even there.

## Screening evidence

42 configurations were screened sequentially at mass 0.1 with 12 RHS, using the same physical operator and full FP64 checks. These were single samples without a timing warmup; the selection used operator work and was confirmed in large-batch finalist/final runs. Screening times are not mixed into the final timing table. Representative configurations:

| Configuration | Aggregate | Vectors/chirality | Filter sweeps | Pre/post | Cycle weight | Coarse dimension | Fine S calls |
|---|---|---:|---:|---:|---:|---:|---:|
| original | 4x4x4x4 | 6 | 8 | 1/1 | 0.65 | 192 | 142 |
| pre_only | 4x4x4x4 | 6 | 8 | 1/0 | 0.65 | 192 | 153 |
| post_only | 4x4x4x4 | 6 | 8 | 0/1 | 0.65 | 192 | 151 |
| filter24 | 4x4x4x4 | 6 | 24 | 0/1 | 0.65 | 192 | 145 |
| time2 | 4x4x4x2 | 6 | 24 | 0/1 | 0.65 | 384 | 105 |
| space2 | 4x4x2x4 | 6 | 24 | 0/1 | 0.65 | 384 | 145 |
| long_x | 8x4x4x2 | 6 | 24 | 0/1 | 0.65 | 192 | 110 |
| basis8 | 4x4x4x4 | 8 | 24 | 1/0 | 0.65 | 256 | 103 |
| basis10 | 4x4x4x4 | 10 | 24 | 1/0 | 0.65 | 320 | 97 |
| basis12_pre | 4x4x4x4 | 12 | 24 | 1/0 | 0.65 | 384 | 83 |
| basis16 | 4x4x4x4 | 16 | 24 | 1/0 | 0.65 | 512 | 83 |
| basis12_filter8 | 4x4x4x4 | 12 | 8 | 1/0 | 0.65 | 384 | 83 |
| basis12_filter96 | 4x4x4x4 | 12 | 96 | 1/0 | 0.65 | 384 | 83 |
| time8_basis24 | 4x4x4x8 | 24 | 24 | 1/0 | 0.65 | 384 | 87 |
| space8_basis24 | 8x4x4x4 | 24 | 24 | 1/0 | 0.65 | 384 | 97 |
| basis12_two | 4x4x4x4 | 12 | 24 | 1/1 | 0.65 | 384 | 82 |
| two_w09 | 4x4x4x4 | 12 | 24 | 1/1 | 0.9 | 384 | 73 |
| two_w10 | 4x4x4x4 | 12 | 24 | 1/1 | 1 | 384 | 73 |
| two_w11 | 4x4x4x4 | 12 | 24 | 1/1 | 1.1 | 384 | 73 |
| two_w13 | 4x4x4x4 | 12 | 24 | 1/1 | 1.3 | 384 | 82 |
| two_w15 | 4x4x4x4 | 12 | 24 | 1/1 | 1.5 | 384 | 149 |
| three_w10 | 4x4x4x4 | 12 | 24 | 2/1 | 1 | 384 | 76 |
| three_pre_w10 | 4x4x4x4 | 12 | 24 | 3/0 | 1 | 384 | 80 |
| two_filter8_w10 | 4x4x4x4 | 12 | 8 | 1/1 | 1 | 384 | 73 |

Removing one smoother alone increased iterations enough to lose the saved Schur work. A 12-vector basis captured the troublesome subspace better; 16 vectors did not reduce these counts further. Shorter time aggregates helped the six-vector basis but did not beat the twelve-vector 4⁴ configuration. Increasing filtering beyond eight sweeps did not improve the selected large-batch convergence. Raising the two-step smoothing weight helped; 1.5 was substantially worse. These are measured choices for this operator, not universal multigrid parameter recommendations.

## Validation, reproduction, and limits

All 96 final measured solves passed every-RHS full CPU FP64 acceptance. Worst final relative residual: 5.31e-11; maximum known-solution relative error: 3.06e-11. The tuned 24-solve small self-test passes, including zero RHS locking. Both masses' actual 384-dimensional tuned and anisotropic hierarchies pass P†P, gamma5 Hermiticity, Galerkin consistency, coarse inverse consistency, CPU/Metal restriction/prolongation, and complete cycle agreement. Pre-only/post-only and larger-basis paths also have saved small/full FP64-checked solves.

```sh
make mg-tune-screen
make mg-tuned-bench
```

The executable's defaults preserve the original settings. To use the selected cycle directly, add `--vectors 12 --omega 1`; see [MULTIGRID.md](MULTIGRID.md) for all options and algorithms. The factory still uses synthetic links and positive masses. The basis is filtered aggregation rather than a production adaptive eigensolver; the coarse matrix uses a CPU FP64 direct inverse, limited to 512 unknowns, with FP32 CPU/GPU applications. Larger lattices, physical ensembles, and near-critical masses need a sparse or iterative coarse solve and separate tuning. No SuperbBLAS or Chroma source is modified.

Raw results: [comparison CSV](results/mg-tuned/comparison.csv), [full FP64 traces](results/mg-tuned/comparison.trace.csv), [final run log](results/mg-tuned-bench.log), [screening](results/mg-screen/). The [metadata note](results/MG_METADATA_NOTE.md) documents correction of a CSV pre-smoothing metadata field from the recorded run arguments; numerical measurements were unaffected.
