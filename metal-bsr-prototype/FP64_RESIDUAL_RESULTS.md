# FP64 residual optimization on the Chroma gauge field

The original residual check spent almost all its time evaluating the full
clover-Wilson stencil. A five-second macOS sample placed the main work in
`CompactClover::apply<double>`. In a nine-pair timing comparison on the
20⁴ field with four resident RHS and eight CPU workers, the original
stencil took a median 193.76 ms and subtraction plus norm accumulation
took 2.27 ms. The replacement took 15.66 ms, approximately 12.5× faster
for the complete check.

The new kernel projects onto the two independent Wilson spin components
before the SU(3) multiplication, then reconstructs the other two components.
It evaluates both chiral clover blocks in double and computes `b - M*x`
and the per-RHS norm during the final stores. It retains two FP64 NEON
lanes for pairs of RHS and handles an odd trailing RHS explicitly.
The working arrays for the stencil and norms are local to each worker.

The optimized and original residual vectors differed by
6.31e-16 relative to `||b||₂` on the profile's dense input and 2.63e-16 on
the saved converged point-source solutions. The latter retained a full
relative residual of 1.261e-11. One- and eight-worker checks cover resident
RHS counts 1, 2, 3, 4, 5, 8, and 12.

The solver continues to recompute the original full equation in FP64 after
every outer correction, for every RHS. Its acceptance threshold is
`||b - M*x||₂ / ||b||₂ <= 1e-10`.

[Baseline sample](results/residual-profile-baseline.sample.txt) ·
[Paired timings](results/residual-profile-comparison.log) ·
[RHS tail checks](results/residual-profile-small.log)

The complete rerun retained the public 20⁴ field, m₀ = −0.375, c_sw = 1,
four resident RHS, eight CPU workers, the original solver settings and
one hierarchy shared by both backends. The 96-RHS result sums the first
24 accepted tiles, and 192 RHS sums all 48.

| RHS | CPU before s | CPU after s | Metal before s | Metal after s |
|---:|---:|---:|---:|---:|
| 96 | 902.98 | 1217.74 | 104.92 | 83.37 |
| 192 | 1850.75 | 2236.73 | 205.72 | 164.97 |

| RHS | CPU FP64 before s | CPU FP64 after s | Metal FP64 before s | Metal FP64 after s |
|---:|---:|---:|---:|---:|
| 96 | 18.16 | 1.83 | 21.44 | 1.71 |
| 192 | 36.72 | 3.45 | 41.90 | 3.30 |

The CPU total increased in this run. At 192 RHS, its measured FP32 inner
work increased from 1813.27 s to 2232.11 s while full FP64 checks fell from
36.72 s to 3.45 s. The inner arithmetic and hierarchy were retained, and
both backends still took four outer corrections, 114–116 fine Krylov steps
and 358–365 fine Schur calls per tile. These sequential runs occurred
under ordinary desktop load. The paired residual microbenchmark above
isolates the cost of the replacement kernel.

All 96 measured tiles passed the full 1e-10 target. The worst final residual
was 1.117e-11. Every final value matches its per-outer trace. A separate
four-point-source solve was written back to disk and independently checked
with Chroma on both backends, reaching worst full residual 1.261e-11.

```sh
make build/chroma_mg build/residual_profile
VECLIB_MAXIMUM_THREADS=1 ./build/residual_profile build/chroma-public20.bridge 9 \
  build/chroma-public20-point4.metal.bin build/chroma-public20-point4.rhs
VECLIB_MAXIMUM_THREADS=1 ./build/chroma_mg \
  --chroma-field build/chroma-public20.bridge --block 5 --tile-rhs 4 \
  --rhs 96,192 --samples 1 --threads 8 --csw 1 --methods mg3 \
  --coarse-iters 8 --bottom-iters 8 --check-hierarchy \
  --csv results/chroma-public20-residual-replay.csv
```

[Measured tiles](results/chroma-public20-residual-optimized.csv) ·
[Full FP64 traces](results/chroma-public20-residual-optimized.csv.trace.csv) ·
[Run log](results/chroma-public20-residual-optimized.log) ·
[Independent Chroma checks](results/residual-optimized-point4.log) ·
[Source and data hashes](results/residual-optimized-artifacts.sha256)
