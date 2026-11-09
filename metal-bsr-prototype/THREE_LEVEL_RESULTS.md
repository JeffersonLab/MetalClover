# Three-level clover-Wilson measurements at 32⁴

Apple M1 Max (24 GPU cores, 32 GB), eight CPU workers, and four resident RHS.
The benchmark uses one paired sweep on the same synthetic gauge field (strength
0.25), with bare mass **m₀ = −0.10**, **c_sw = 1**, and antiperiodic time.
The mass is on the positive side of the low-gap region found by the
[spectral scan](results/compact-critical32.csv). The scan did not calibrate
a physical critical mass.

## Measured solve time

Times are minutes. Each prefix sums measured four-RHS tiles from offset zero.
CPU smoother-only FGMRES reached 48 RHS before the run was paused, so 48 RHS
is the largest matched comparison across all six cases.

| Solver | CPU 48 | Metal 48 | CPU 96 | Metal 96 | CPU 192 | Metal 192 |
|---|---:|---:|---:|---:|---:|---:|
| Three-level MG | 16.17 | 3.60 | 31.89 | 6.21 | 65.26 | 11.26 |
| BiCGStab | 81.64 | 9.67 | 167.50 | 19.08 | 362.37 | 35.53 |
| Smoother-only FGMRES | 241.82 | 13.59 | — | 29.85 | — | 66.35 |

At 192 RHS, Metal MG is **5.80× faster than CPU MG**, **3.16× faster than
Metal BiCGStab**, and **5.89× faster than Metal smoother-only FGMRES**.
Of the 11.26 Metal MG minutes, 6.40 were FP32 inner work and 4.73 were
full CPU FP64 residual checks. Peak recorded Metal allocation for the
four-RHS MG case was 8.48 GiB.

Solve time covers even/odd reduction and reconstruction, Krylov work, FP64
updates, and the full residual checks. Setup, warmup, manufactured RHS
generation, and the independent known-solution check are outside solve time.
The run paused for a reboot and resumed from accepted tiles with a fresh
warmup. Wall times include about 17 minutes between lid close and full wake
in the CPU BiCGStab 192-RHS total and about 42 minutes in the CPU smoother
48-RHS total. [Checkpoint and power log notes](RESUME_THREE_LEVEL.md).

## Convergence

All **252 accepted tiles** (1,008 RHS executions) passed the full CPU FP64
relative residual target of **1e-10** on the original reconstructed equation.
The worst accepted residual was **9.90e-12**. The audit matched every accepted
tile to its final per-outer trace. Three-level MG needed four outer
corrections, 29 fine Krylov steps and 91 fine Schur calls per tile on both
processors.

| Solver | CPU fine Schur calls per tile, first 48 RHS | Metal fine Schur calls per tile, first 48 RHS |
|---|---:|---:|
| Three-level MG | 91 | 91 |
| BiCGStab | 534–782 | 644–1,848 |
| Smoother-only FGMRES | 1,311–1,673 | 1,224–1,688 |

An earlier four-RHS Metal pilot at m₀ = −0.12 passed with three-level MG in
94 fine Schur calls and a full residual of 8.4e-13. Plain Metal BiCGStab hit
its 2,048-step inner cap at that mass before producing an accepted full
solution. The timing comparison uses m₀ = −0.10 for every case.

[Methodology and reproduction commands](THREE_LEVEL.md) ·
[Raw measurements](results/threelevel32-compare.csv) ·
[Per-outer FP64 traces](results/threelevel32-compare.csv.trace.csv)
