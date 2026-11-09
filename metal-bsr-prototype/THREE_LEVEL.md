# Three-level clover-Wilson prototype

This extension stays inside the standalone prototype. Chroma and SuperbBLAS
are unchanged. It retains `c_sw=1`, the QDP gamma basis with diagonal gamma5,
two local Hermitian 6x6 clover blocks per site, and the odd Schur equation

```
S_o = A_oo - H_oe A_ee^-1 H_eo,
M   = (m0+4)I - Dslash/2 + C.
```

The gauge-strength parameter remains 0.25, with antiperiodic time. The gauge
generator is the same as the earlier experiments. These are synthetic fields.

## Hierarchy

The first interpolation has 12 vectors per chirality on 4^4 aggregates of
the odd fine grid. At 32^4 this produces 4,096 coarse sites with 24 components
each: **98,304 unknowns**. The second interpolation uses 12 vectors per chirality
on 2^4 aggregates of that 8^4 coarse lattice, producing **6,144 unknowns**.
Both bases start from independent constant-plus-noise test vectors, undergo
eight filtering sweeps at weight 0.65, and use two-pass local orthogonalization
within each chirality. Neither basis uses the manufactured solution or RHS.

Both coarse matrices are sparse exact Galerkin projections. The first is
assembled through local even-site contributions without a global dense
interpolation matrix or a volume-sized intermediate for every aggregate pair.
The second is assembled by sparse block contractions. Their FP32 coefficient
storage is 594 MiB and 37.125 MiB respectively. Global dense inverses are not
constructed; only the local 6x6 and 24x24 diagonal blocks are inverted.

Fine FGMRES uses restart 8. Its three-level preconditioner has one pre/post
clover-block Jacobi step at weight 1. The first coarse equation uses FGMRES(4),
relative target 0.1 and an eight-iteration cap. Its preconditioner has one
pre/post block-Jacobi step and the second coarse correction. The bottom solve
is block-Jacobi-preconditioned BiCGStab with target 0.1 and cap 8. Its target uses the left-preconditioned residual norm; the actual
unpreconditioned bottom residual is a separate diagnostic. A capped
coarse solve can return an inexact correction; flexible GMRES stores that
correction. Metal's bottom solver uses per-RHS masks within its fixed dispatch
budget. CPU stops submitting products once it converges. These costs are
included in the measurements.

The CPU fine BiCGStab reduction uses FP32 sums over blocks of 256 components,
with a final FP64 sum of those partials, cast back to FP32 scalar values.
This avoids very long unbroken FP32 dot-product accumulations at 32^4.
Metal reductions remain FP32 on the GPU.

The baselines remain BiCGStab on the same odd Schur equation and FGMRES(8)
with two clover-block Jacobi steps. Fine inner target is 0.002, with a
2,048-iteration safety cap. Every accepted full solution must meet 1e-10
on every RHS in the original reconstructed equation.

## Compact operator and memory

`compact_clover.hpp` stores forward links, neighbor indices and local clover
blocks instead of three expanded copies of the 17-term matrix and its even/odd
submatrices. CPU and Metal apply the same Wilson projectors. Small-lattice
tests compare this FP64 operator and its Schur complement against the original
encoded implementation, independently of the new sparse assembly.

The outer solve retains three full FP64 fields: b, x, and the residual.
The full residual buffer is overwritten by `M64*x` after each correction.
The known solution is regenerated from the deterministic component/RHS hash
for its independent error check, which follows acceptance outside solve time.
It is not kept in another full field. CPU correction buffers and Metal upload/
reconstruction buffers are reused after their previous contents are no longer
needed. Metal even/odd preparation and reconstruction remain FP32 GPU work.

Full FP64 gauge links and clover blocks remain resident for every residual
check. FP64 interpolation/coarse setup copies are released after hierarchy
validation; the working hierarchy is FP32 on both processors. Small QR
problems use CPU FP64. No FP64 GPU arithmetic is assumed.

An eight-RHS trial allocated about 16 GiB in Metal alone, in addition to
the CPU FP64 fields and operator data. The full study uses a smaller tile
to leave headroom for those fields and ordinary desktop activity.

The intended resident tile is **four RHS**. RHS 0–95 and 0–191 use 24 and
48 actual sequential tiles. CPU and Metal receive identical global RHS IDs.
The 96-RHS result is the measured prefix of the 192-RHS sweep, not an
extrapolation. Every tile starts from zero. Checkpoints record accepted tiles,
and `--resume` skips them after checking numerical configuration and acceptance.

## Mass selection

The separate compact probe minimizes Rayleigh values of `Q^2=M†M`, with
`Q=gamma5*M`, in the span of the current vector, its residual and the previous
search direction. It keeps a fixed number of full vectors rather than storing
a growing Lanczos basis. It recomputes the actual full-lattice FP64 eigen-
residual and checks free-field estimates against the analytic dispersion.

The 32^4 scan covers bare masses -0.10, -0.12 and -0.14, with an independent
longer start at -0.12. Rayleigh gap estimates are around 0.09. These runs do
not attain their requested eigen-residual tolerance. Their upper estimates
and stationary-mass diagnostics locate a low-gap region but do not certify
the smallest singular value or calibrate a physical PCAC critical mass.

The -0.12 pilot passed with three-level MG (30 Krylov steps, 94 fine Schur
applications, full residual about 8.4e-13 for four RHS). Plain Metal BiCGStab
did not meet its inner target within 2,048 iterations on that pilot. The
comparison therefore tests **m0=-0.10**, toward the positive side of the same
low-gap region, while retaining all solver and full-residual targets.

## Measurement protocol

The full study uses one paired sweep per solver and processor, with eight CPU
workers and single-threaded setup BLAS. One complete resident-tile warmup
precedes each measured method. Metal methods run sequentially, followed by
CPU methods. Large method workspaces are released between methods. Common
operator and hierarchy setup is retained across the study.

Solve time includes initialization, FP32 elimination/reconstruction, all fine
and coarse work, vector operations, CPU small QR, GPU submissions/waits,
FP64 solution updates and full CPU FP64 residual checks. It excludes setup,
manufactured RHS construction, warmup and the independent known-solution
error check. Per-outer progress logging is inside total solve time. Ordinary
local system activity is present; there is no thermal/scheduling isolation.
The pilot and functional-check times are not final benchmark samples.

## Commands

```sh
make build/metal_three_level build/compact_probe
VECLIB_MAXIMUM_THREADS=1 ./build/metal_three_level --self-test \
  --mass -.10 --csv results/threelevel-self-test-final.csv
VECLIB_MAXIMUM_THREADS=1 ./build/metal_three_level --size 32 \
  --mass -.10 --validate-only
VECLIB_MAXIMUM_THREADS=1 caffeinate -i ./build/metal_three_level \
  --size 32 --mass -.10 --csw 1 --tile-rhs 4 --rhs 96,192 --samples 1 \
  --threads 8 --coarse-iters 8 --bottom-iters 8 \
  --csv results/threelevel32-compare.csv
```

The same benchmark command with `--resume` appends only missing tiles.
Accepted CSV rows and the per-outer full FP64 trace are flushed for each tile.
`python3 threelevel_progress.py` prints a read-only progress snapshot and
completed measured totals for the 96/192-RHS prefixes.

Primary files: `compact_clover.hpp`, `compact_clover.metal`,
`compact_solve.hpp`, `three_level.hpp`, `three_level.mm`, `compact_probe.mm`.
The original FGMRES/coarse solver helpers now accept generic operator types;
their existing two-level aliases remain available and pass regression checks.
