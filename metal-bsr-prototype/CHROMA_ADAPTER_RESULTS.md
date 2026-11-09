# Chroma adapter benchmark on a 20⁴ gauge configuration

The benchmark used the public SIMULATeQCD NERSC archive
`l20t20b06498a_nersc.302500`. Chroma verified its archive checksum and
plaquette. The header gives an average plaquette of 0.6381995717.
The valence operator was isotropic clover-Wilson at **m₀ = −0.375**,
**c_sw = 1**, **u₀ = 1**, with periodic spatial and antiperiodic temporal
fermion phases. The mass is a benchmark parameter on this ensemble.
[Public gauge file](https://github.com/LatticeQCD/SIMULATeQCD/tree/main/test_conf).

## Chroma operator match

Chroma 3.44.0 loaded the NERSC links and computed the FP64 reference actions
from its QDP++ clover term and Wilson dslash. The standalone code built its
own two chiral 6×6 clover blocks and applied the odd Schur complement
`A_oo - H_oe A_ee^-1 H_eo`. Relative 2-norm differences on a dense Z4
probe were **1.51e-16** for the full operator and **1.60e-16** for the odd
Schur complement. The two Galerkin operators matched direct FP64
projections within 2.51e-15. CPU and Metal transfer and coarse sparse
products agreed within 8.90e-7.

## Measured time on M1 Max

One paired sweep with eight CPU workers, four RHS resident, and independent
deterministic Z4 sources. Both backends used the same hierarchy, with 5⁴
fine aggregates followed by 2⁴ coarse aggregates and 12 vectors per
chirality at each interpolation. Odd fine components number 960,000,
then 6,144 and 384 at the two coarse levels. The 96-RHS result is the first
24 measured tiles and the 192-RHS result sums all 48.

| RHS | CPU solve s | Metal solve s | CPU / Metal | Full CPU FP64 checks inside Metal s |
|---:|---:|---:|---:|---:|
| 96 | 902.98 | 104.92 | 8.61× | 21.44 |
| 192 | 1850.75 | 205.72 | 9.00× | 41.90 |

Tile time includes Z4 RHS preparation, even/odd reduction and reconstruction,
FP32 Krylov work, FP64 updates and full FP64 checks. It excludes archive
export, one-time initialization, the reference checks, one warmup tile per
backend, and file output. The run measured 0.035 s for bridge link import,
0.190 s for compact operator setup, and 6.057 s for the one shared hierarchy
build. Largest recorded Metal allocation during the Metal solves was 1.38 GiB.
This was one sequential run on a laptop under ordinary desktop load.

Every RHS passed `||b-Mx||₂ / ||b||₂ ≤ 1e-10` on the original full equation.
The worst residual across both sweeps was **1.117e-11**. On both backends
each tile needed four outer corrections, 114–116 fine Krylov steps and
358–365 fine Schur applications.

## File round trip

On the same 20⁴ field, a separate four-RHS tile loaded point sources from
disk and wrote the CPU and Metal solutions back in the bridge's spinor order.
Chroma independently re-applied its FP64 full operator to every returned
solution. The worst residual was **1.261e-11** on each backend.
[File round-trip log](results/chroma-public20-point4.log).

## Reproduce

The archive SHA-256 is
`26a261aa730ef6e4f1ca06b3e7d1cd76bbecf622222490c73db79f8f39356bdd`.

```sh
mkdir -p build
curl -fL --output build/l20t20b06498a_nersc.302500 \
  https://media.githubusercontent.com/media/LatticeQCD/SIMULATeQCD/main/test_conf/l20t20b06498a_nersc.302500
OMP_NUM_THREADS=1 ./build/chroma_export build/l20t20b06498a_nersc.302500 \
  build/chroma-public20.bridge -.375
VECLIB_MAXIMUM_THREADS=1 ./build/chroma_mg \
  --chroma-field build/chroma-public20.bridge --block 5 --tile-rhs 4 \
  --rhs 96,192 --threads 8 --csw 1 --methods mg3 --coarse-iters 8 \
  --bottom-iters 8 --check-hierarchy --csv results/chroma-public20-replay.csv
```

[Build and file layout](CHROMA_ADAPTER.md) ·
[Measured tiles](results/chroma-public20.csv) ·
[FP64 residual traces](results/chroma-public20.csv.trace.csv) ·
[Run log](results/chroma-public20.log) ·
[Artifact hashes](results/chroma-public20-artifacts.sha256)
