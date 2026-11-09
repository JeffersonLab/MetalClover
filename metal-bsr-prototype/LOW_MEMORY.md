# Fitting the real 32³×256 clover operator

The low-memory driver retains the exact full FP64 action for defect checks.
GPU local blocks use `A = L D L†`, with unit-diagonal lower `L`. Each
chirality stores six real diagonal entries and fifteen complex lower entries
(36 real values), rather than 36 complex matrix entries. The same odd factor
is used for `A_oo*x` and triangular `A_oo^-1` solves in the smoother. The
even factor is used for Schur elimination and RHS reconstruction. There is
no separate GPU odd inverse array. Positivity of every local LDL pivot is
checked during FP64 setup for this action.

The CPU gauge is mapped from the native FP64 stout bridge. CPU clover blocks
are regenerated from those original FP64 links at each site, with a small
per-thread cache; they do not come from the GPU factors. Only the full
solution `x` is persisted in FP64. The entire full-system residual is
evaluated in FP64 in eight-time-slice slabs, including its norm, before
casting each slab into the GPU's existing FP32 correction input. Its FP64
values need not be retained afterward. Source values are regenerated from
the same deterministic complex Z4 hash used by the native verifier.

Dirty solution slabs are synchronously written back before their pages are
dropped. Gauge pages are clean, file-backed, and dropped between slabs.
The eight-time-slice double residual scratch is anonymous and discarded
after use. Full FP64 residual checks and their clover recomputation cost are
included in solve timings. Each field must reach `||b-Mx||/||b|| <= 1e-12`.
An independent executable then checks every saved field with original
Chroma `UnprecCloverLinOp` in FP64.

`CHRMG003` adds four padding bytes to the existing 76-byte bridge header,
aligning FP64 gauge/probe data at byte 80 for mapping. Both bridge versions
remain supported by the ordinary streaming reader. Native export selects
the aligned format with `METALCLOVER_ALIGNED_BRIDGE=1`. Physical gauge links
remain unweighted FP64; anisotropy and antiperiodic time are applied once
by the respective operator implementations.

```sh
make -f Makefile.production build/low_memory_solve build/low_memory_verify \
  CHROMA_CONFIG=/path/to/chroma-config CHROMA_SOURCE=/path/to/chroma
VECLIB_MAXIMUM_THREADS=1 /private/tmp/low_memory_solve \
  BRIDGE_PREFIX OUTPUT_PREFIX mg 2 4,4,4,8 8 128
# Use bicg instead of mg for the separate GPU BiCGStab comparison.
OMP_NUM_THREADS=8 /private/tmp/low_memory_verify \
  ACTION_XML BRIDGE_PREFIX OUTPUT_PREFIX metal-mg 2 CHECKS_CSV
```

Stage executables outside FileProvider before running. Use a local scratch
prefix for the large mapped bridge and fields. The two solver processes
run sequentially; one MG hierarchy is constructed and reused across its
RHS loop. The GPU fine FGMRES restart is 8; both transfers have eight basis
vectors per chirality, with fine aggregates 4×4×4×8 and second aggregates
2×2×2×2. Fine filtering uses 128 Metal sweeps, second-level filtering 16
CPU FP64 sweeps. One pre/post block-Jacobi step, coarse FGMRES restart 8 /
limit 32 / target 0.1, and bottom BiCGStab limit 4 / target 0.1 follow the
24³ tuning, with a shorter fine restart to fit the larger lattice.

The real header action is mass -0.0860, xi0 4.3, nu 1.265, clover R
1.58932722549812 / T 0.902783591772098, two spatial stout steps at rho 0.14,
and fermion boundary `(1,1,1,-1)`. The source tags are 20261107 and 20261108.
