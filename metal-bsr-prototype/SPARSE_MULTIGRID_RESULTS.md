# Sparse coarse multigrid: 16⁴ CPU versus Metal

Local measurements on Apple M1 Max (24 GPU cores, 32 GB, macOS 26.7.1). The standalone clover-Wilson operator keeps c_sw=1, QDP gamma conventions, gauge strength 0.25, positive bare masses 0.1 and 0.5, and antiperiodic time. The comparison uses the same sparse two-level solver on CPU and Metal, with complete CPU FP64 residual checks on the reconstructed full equation for every RHS.

## Complete solve time

Medians of three measured logical batches, from zero. **96 RHS are resident simultaneously. 192 RHS run as two sequential 96-RHS tiles**, sharing the hierarchy/operator/inner workspace. CPU and Metal workspaces are allocated sequentially to fit the 32 GB Mac. One complete resident-tile warmup precedes each backend's timings for each mass. Backend order reverses between masses. CPU uses eight GCD workers and single-threaded Accelerate for setup. The computer had ordinary local activity, without thermal or scheduler isolation.

| Mass | Logical RHS | Resident RHS | CPU total s | Metal total s | CPU / Metal | CPU RHS/s | Metal RHS/s |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 96 | 33.538 | 7.838 | 4.28× | 2.86 | 12.25 |
| 0.1 | 192 | 96 | 75.889 | 15.767 | 4.81× | 2.53 | 12.18 |
| 0.5 | 96 | 96 | 26.364 | 6.478 | 4.07× | 3.64 | 14.82 |
| 0.5 | 192 | 96 | 48.270 | 13.068 | 3.69× | 3.98 | 14.69 |

Solve totals include resets, FP32 RHS reduction and even-site reconstruction, fine and coarse sparse products, smoothing, FGMRES/BiCGStab vector operations and reductions, CPU small QR, GPU submissions/waits, FP64 accumulation, and full FP64 residual checks. Synthetic manufactured RHS generation is outside solve time and has its own raw column; supplied RHS would not require that operator application. Hierarchy, problem, and inner workspace setup are reported separately. FP64 outer buffers are allocated once per logical case outside the solve timer; their allocation time was not collected. Reset and data movement inside solve are timed. No 192-RHS result represents 192 simultaneously resident Krylov vectors.

| Mass | RHS | CPU min–max s | Metal min–max s | Metal CPU-FP64 residual s | Metal FP32 inner s |
|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 32.748–39.710 | 7.756–7.850 | 2.317 | 5.400 |
| 0.1 | 192 | 64.681–76.187 | 15.666–15.884 | 4.760 | 10.766 |
| 0.5 | 96 | 24.946–29.492 | 6.464–6.574 | 2.322 | 4.059 |
| 0.5 | 192 | 48.039–48.696 | 12.923–14.474 | 4.613 | 8.209 |

## Convergence and work

The fine solver is FGMRES(8) with an inner target of 0.002. The basis has twelve vectors per chirality, 4⁴ aggregates, eight FP64 filtering steps at weight 0.65, and one pre/post clover-block Jacobi step at weight 1.0. Fine S counts include smoothing, Arnoldi applications, and actual FP32 residual refreshes. Counts are summed over resident tiles for the logical batch; a call always applies to 96 resident RHS.

| Mass | RHS | Backend | Outer corrections/tile | Fine Krylov steps/tile | Total fine S calls | Total coarse SpMV dispatches | Total coarse iterations | Worst full FP64 residual |
|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | cpu_sparse_mg | 4 | 24 | 76 | 135 | 60 | 9.59e-13 |
| 0.1 | 96 | metal_sparse_mg | 4 | 24 | 76 | 216 | 60 | 9.59e-13 |
| 0.1 | 192 | cpu_sparse_mg | 4 | 24 | 152 | 269 | 119 | 9.68e-13 |
| 0.1 | 192 | metal_sparse_mg | 4 | 24 | 152 | 432 | 119 | 9.68e-13 |
| 0.5 | 96 | cpu_sparse_mg | 4 | 19 | 61 | 59 | 20 | 2.38e-13 |
| 0.5 | 96 | metal_sparse_mg | 4 | 19 | 61 | 171 | 20 | 2.38e-13 |
| 0.5 | 192 | cpu_sparse_mg | 4 | 19 | 122 | 118 | 40 | 2.39e-13 |
| 0.5 | 192 | metal_sparse_mg | 4 | 19 | 122 | 342 | 40 | 2.39e-13 |

The coarse solver is left block-Jacobi-preconditioned complex BiCGStab, starting from zero, with target relative preconditioned residual 0.1 and a four-iteration cap. The returned coarse correction may be inexact; FGMRES stores Z to accommodate that variable preconditioner. The unpreconditioned FP32 coarse residual is recomputed after every correction. It is a diagnostic and can exceed 0.1. CPU stops submitting sparse products after convergence; Metal encodes the cap with per-RHS masking. GPU SpMV counts include masked/padded dispatches and both backends include one true-coarse-residual product per correction. This difference is included in the runtime comparison.

| Mass | RHS | Worst CPU coarse true residual | Worst Metal coarse true residual |
|---:|---:|---:|---:|
| 0.1 | 96 | 0.1268 | 0.1267 |
| 0.1 | 192 | 0.1966 | 0.1961 |
| 0.5 | 96 | 0.0993 | 0.0993 |
| 0.5 | 192 | 0.0993 | 0.0993 |

Example full FP64 convergence history, mass 0.1, 192 RHS, trial 0:

| Global RHS offset | Outer correction | CPU full residual | Metal full residual |
|---:|---:|---:|---:|
| 0 | 1 | 0.000513 | 0.000513 |
| 0 | 2 | 8.66e-07 | 8.66e-07 |
| 0 | 3 | 1.7e-09 | 1.7e-09 |
| 0 | 4 | 9.59e-13 | 9.59e-13 |
| 96 | 1 | 0.000514 | 0.000514 |
| 96 | 2 | 8.65e-07 | 8.65e-07 |
| 96 | 3 | 1.68e-09 | 1.68e-09 |
| 96 | 4 | 9.68e-13 | 9.68e-13 |

## Sparse storage, setup, and memory

At 16⁴ the retained odd space has 393,216 complex components. 256 aggregates × 24 components give 6,144 coarse unknowns. The exact Galerkin operator has 8,448 dense 24×24 blocks (33 per aggregate), stored as block CSR. Its FP32 coefficients occupy **37.125 MiB**, compared with 288 MiB for a full FP32 coarse matrix. An FP64 sparse copy, sparse indices, basis, and local diagonal inverses are additional. Only local 24×24 blocks are inverted; no global dense inverse, factorization, or fine-by-all-coarse prolongator array is formed. The latter would take 36 GiB per FP64 array at this geometry, before all other workspaces.

| Mass | RHS | Hierarchy setup s | CPU inner workspace setup s | Metal inner workspace setup s | Problem setup s | Process peak RSS GiB |
|---:|---:|---:|---:|---:|---:|---:|
| 0.1 | 96 | 2.757 | 1.101 | 0.170 | 2.645 | 16.34 |
| 0.1 | 192 | 2.757 | 1.101 | 0.170 | 2.645 | 16.34 |
| 0.5 | 96 | 3.961 | 0.807 | 0.063 | 2.490 | 18.12 |
| 0.5 | 192 | 3.961 | 0.807 | 0.063 | 2.490 | 18.12 |

Hierarchy setup includes FP64 test-vector filtering, local QR, sparse assembly, local block inversion, and temporary setup operator construction. One mass-specific hierarchy is reused across both processors, 96/192 RHS, and all trials. Process RSS is the cumulative OS high-water mark within this run; it is not a per-method isolated measurement or total system memory usage. Unrestricted simultaneous 192-RHS fine Krylov stacks would impose a much larger memory footprint. The user-visible 192-RHS batch therefore uses two bounded-memory tiles. With many RHS, hierarchy setup is paid once per gauge/mass/geometry, not once per tile. The setup table is not a complete cold-start total because FP64 outer-buffer allocation was not collected.

## Validation and reproduction

All 24 measured logical batches passed the full FP64 criterion for every RHS. They contain 36 resident-tile solves and 3456 RHS executions across repeated trials. Worst final true relative residual was 9.68e-13; maximum known-solution relative error was 2.62e-12. The trace records each tile's global RHS offset and every outer true residual. A known solution is used only for validation and RHS generation; hierarchy setup is independent of it.

Before timings, the 40-solve clover regression, 24-solve dense multigrid regression, and eight-solve sparse self-test passed. Sparse self-tests compare with independent dense column probing, coarse gamma5 Hermiticity, CPU/Metal sparse products, and strict coarse solves, including zero RHS locking. Actual 8⁴ and 16⁴ hierarchies passed FP64 Galerkin/coarse-solver checks at both masses. At 16⁴ the Galerkin random-vector identity agrees at about 2e-15 relative error.

```sh
make build/metal_sparse_multigrid build/metal_multigrid build/metal_clover
python3 run_sparse_study.py
python3 sparse_summary.py
```

The main timing command alone is:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/metal_sparse_multigrid \
  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass .1,.5 --csw 1 \
  --vectors 12 --aggregate 4 --restart 8 --coarse-iters 4 \
  --coarse-tol .1 --samples 3 --threads 8 --csv results/sparse-mg16.csv
```

This establishes CPU/Metal convergence and timing for the sparse iterative coarse solver at the larger volume. The gauge field remains synthetic and the masses positive; production ensembles, critical-mass behavior, distributed execution, and additional multigrid levels need separate validation. See [SPARSE_MULTIGRID.md](SPARSE_MULTIGRID.md) for the assembly equations and implementation details.

Raw data: [results/sparse-mg16.csv](results/sparse-mg16.csv), [full FP64 traces](results/sparse-mg16.csv.trace.csv), [benchmark log](results/sparse-mg16.log), [regression/study log](results/sparse-study.log).
