# Fitting the real 32³×256 clover operator

## Completed real-field comparison

The real **32³×256, mass -0.0860** field now fits and solves on the **32 GiB
Apple M1 Max** (macOS 26.7.1). Two deterministic dense Z4 RHS were solved from
zero, one resident at a time. One MG hierarchy was built and reused for both.
All four saved GPU fields passed the independent original Chroma FP64
full-system check at `1e-12`.

| Method | Mean solve (s/RHS) | Fine Schur applications/RHS | GPU allocation at ready (GiB) | Maximum original Chroma FP64 residual |
|---|---:|---:|---:|---:|
| Metal BiCGStab | 2211.211 | 12714.0 | 9.285 | 2.5931e-14 |
| Tuned Metal MG / FGMRES | **916.269** | **2743.5** | **22.700** | **2.4821e-14** |

MG is **2.41× faster** on solve time, with **78.4% less fine Schur
work**. BiCGStab took **2192.742 / 2229.680 s**; MG took **893.415 / 939.123 s**.
BiCGStab used **6263 / 6451** iterations and **12526 / 12902** fine Schur
applications. MG used **880 / 874** FGMRES iterations and **2753 / 2734**
fine applications. All four solves used five defect corrections; their
final full FP64 residuals are at comparable accuracy.

The single MG hierarchy cost **269.226 s** to construct and upload, including
workspace setup. Amortized over the two tested RHS, MG took
**1050.882 s/RHS**, still **2.10× faster**. Based on these two-source
measurements, arithmetic setup amortization is
`916.269253 + 269.226134/N` seconds/RHS. Setup adds 2.804 s/RHS at N=96 and
1.402 s/RHS at N=192. These are projections; this comparison actually
measured two RHS, with one resident at a time. Setup filtering applied
2048 additional fine Schur operators, counted in setup rather than solves.

The measured **22.700 GiB** GPU allocation is below the device's
**24.960 GiB** recommended working set. The CPU full residual scratch is
**48 MiB**, plus halo/gauge mapping pages and geometry data. The complete
FP64 solution is a 1.5 GiB file-backed field. Full CPU clover, source, and
residual arrays are absent. Clean mapping pages are dropped between slabs;
only changed solution slabs require synchronous writeback. The reported
Metal allocation is authoritative: the primary-buffer accounting predicted
22.325 GiB, with a further 0.375 GiB in measured allocation not disaggregated
between allocator/runtime/setup-retained costs. Original native verification
ran separately, with peak RSS 17.944 GiB for BiCGStab fields and 18.742 GiB
for MG fields. No inversion and native verification overlapped.

The GPU local clover storage fell from **6.75 GiB** (odd A, even inverse,
odd inverse as full complex matrices) to **2.25 GiB** (one LDL factor per
parity, with the odd factor shared by Schur multiplication and smoothing),
saving **4.50 GiB**. A shorter fine restart also reduces Krylov storage.
The large CPU FP64 storage reduction comes from exact recomputation and
streaming rather than promotion of GPU coefficients.

Every inline full residual is computed in CPU FP64 from the original native
FP64 stout gauge, including newly regenerated FP64 clover blocks. Its
unnormalized dense full-operator probe agrees with exported original Chroma
to **1.7086e-16** relative on the entire lattice. The separate original
`UnprecCloverLinOp` probe agrees exactly with that export. Final native MG
residuals are **2.4630e-14 / 2.4821e-14**; BiCGStab residuals are
**2.5931e-14 / 1.9624e-14**. No original CPU inversion was run on 32³×256;
the CPU reference here is original Chroma's full operator evaluation of
every final field.

MG's inline FP64 checks cost **212.391 s** over both RHS
(**11.59%** of solve time); BiCGStab's cost
**126.028 s**. Both costs are inside the table's solve times, as are
FP64 updates, solution-slab writeback, and Metal RHS reduction/reconstruction.
Input import, common operator setup/probe, source generation, shader
compilation, and subsequent independent native verification are outside
solve timing. The hierarchy cost is separately included above. Fine Schur
counts include fine smoothing, outer products and FP32 residual refreshes;
RHS/reconstruction hopping and CPU FP64 checks are inside timing but outside
that particular count.

The first coarse FGMRES solve averaged **3.91 iterations** per call and never
hit its cap of 32. Coarse solves remain inexact preconditioner applications;
the bottom target 0.1 is for its left-preconditioned system. The separately
recorded worst bottom unpreconditioned FP32 ratio was 2.121. That value does
not set acceptance of the physical solution: every saved field must pass
the full CPU FP64 criterion and then original Chroma verification.

The completed final-path preflight used two real 12³×96 RHS for each solver;
all four independent native residuals were below 1.84e-14. The stronger
256/64-sweep, bottom-cap-8 MG configuration then passed an additional small
full solve with native residual 1.1304e-14 before the large MG run.

[Structured results](results/performance32-fit/summary.json),
[BiCGStab timings](results/performance32-fit/test.bicg.csv),
[MG timings](results/performance32-fit/test.mg.csv),
[BiCGStab native checks](results/performance32-fit/test.bicg-native.csv),
[MG native checks](results/performance32-fit/test.mg-native.csv),
[measured memory](results/performance32-fit/mg-measured-memory.json), and
[accepted field/bridge hashes](results/performance32-fit/accepted-files.sha256)
record the evidence. Full traces, resource samples, native logs, exact header
action, runtime file hashes, source hashes, and executable hashes are in
that result directory. Large binary fields remain generated local scratch
artifacts, excluded from Git.

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
  BRIDGE_PREFIX OUTPUT_PREFIX mg 2 4,4,4,4 8 256 64 8
# Use bicg instead of mg for the separate GPU BiCGStab comparison.
OMP_NUM_THREADS=8 /private/tmp/low_memory_verify \
  ACTION_XML BRIDGE_PREFIX OUTPUT_PREFIX metal-mg 2 CHECKS_CSV
```

Stage executables outside FileProvider before running. Use a local scratch
prefix for the large mapped bridge and fields. The two solver processes
run sequentially; one MG hierarchy is constructed and reused across its
RHS loop. The GPU fine FGMRES restart is 8; both transfers have eight basis
vectors per chirality, with fine aggregates 4×4×4×4 and second aggregates
2×2×2×2. Fine filtering uses 256 Metal sweeps, second-level filtering 64
CPU FP64 sweeps. One pre/post block-Jacobi step, coarse FGMRES restart 8 /
limit 32 / target 0.1, and bottom BiCGStab limit 8 / target 0.1 follow the
24³ tuning, with a shorter fine restart, smaller temporal aggregates, and
stronger basis filtering for this larger, lighter configuration. The optional
last two CLI arguments set second-level filter sweeps and the bottom limit;
omitting them preserves the original 16/4 defaults.

The real header action is mass -0.0860, xi0 4.3, nu 1.265, clover R
1.58932722549812 / T 0.902783591772098, two spatial stout steps at rho 0.14,
and fermion boundary `(1,1,1,-1)`. The source tags are 20261107 and 20261108.
