# Homogeneous near-null solves and chiral SVD compression

`low_memory_solve` now has an optional `--basis near-null-svd` builder inspired
by Eloy Romero's `getMGProlongator` in `mgproton.cc`. The original filtered
builder remains the default and is selected explicitly by `--basis filtered`.
The fine action, Schur operator, local LDL factors, physical RHS, outer
refinement, and original Chroma verifier are shared by both modes.

## Real 24³×256 comparison

Two deterministic dense Z4 RHS were solved on the same mass **-0.0856**
configuration and header action as [the earlier tuning study](MG_TUNING_REAL.md).
One RHS was resident at a time. Each MG method constructed **one hierarchy**
and reused it for both sources. All six saved fields passed the original
Chroma full-system FP64 check at the requested **1e-12** tolerance.

| Method | Mean solve (s/RHS) | One hierarchy (s) | Fine Schur applications/RHS | Maximum original Chroma FP64 residual |
|---|---:|---:|---:|---:|
| **Metal MG, near-null/SVD** | **85.561** | 426.623 | **420.5** | 1.5882e-14 |
| Metal MG, filtered control | 260.926 | 65.207 | 2318.5 | 2.3846e-14 |
| Fresh Metal BiCGStab | 671.971 | — | 9072.0 | 2.2605e-14 |

The new basis gives **7.85× faster warm solves than Metal BiCGStab** and
**3.05× faster warm solves than filtered MG**. Fine Schur work falls by
**95.36%** and **81.86%**, respectively. Individual solve times were
84.900 / 86.221 s for the new basis, 263.919 / 257.932 s for filtered MG, and
662.159 / 681.782 s for BiCGStab. Their Krylov iteration counts were
134 / 137, 751 / 752, and 4453 / 4619. Every solve used five defect corrections.

The setup tradeoff matters. Including the one hierarchy over the two measured
RHS gives **298.872 s/RHS** for near-null/SVD MG, **293.529 s/RHS** for filtered
MG, and **671.971 s/RHS** for BiCGStab. Thus the new method is already **2.25×
faster than BiCGStab including setup**, while filtered MG is slightly faster
for just two RHS. The measured means put the crossover against filtered MG
at 2.06 RHS, or **three RHS** after rounding up.

For many sequential RHS, the new method's arithmetic amortization is
`85.5605568745 + 426.62347625/N` seconds/RHS. This gives **90.005 s/RHS at
N=96** and **87.783 s/RHS at N=192**. Only two RHS were measured: these larger
counts project the same mean single-RHS cost and one hierarchy, without a
claim of measured batched throughput.

The 24 candidates were compressed to eight vectors per chirality on
3×3×3×8 aggregates. SVD retained **97.91%** of candidate energy globally and
at least **92.61%** in every aggregate/chirality. The smallest retained
singular-value ratio was 0.0326, and the worst local FP64 orthogonality error
was **6.22e-15**. All 24 setup solves hit their 128-step cap; their measured
`||S*v||/||S*y||` ratios were 0.0192–0.0207, relative to the warmed starting
vectors. They are inexact setup candidates, separate from physical acceptance.
Candidate generation dominates setup; the compression stage took **5.473 s**.
Setup used 4152 fine Schur applications, versus 1024 for the filtered control.

The full operator probe matched the original Chroma export to **1.342e-16**;
the independent native operator probe matched to **1.739e-16**. CPU FP64 inline
residual checks consumed **38.795 s** across the two new-MG solves, or **22.7%**
of their solve time. Those checks, even/odd reduction/reconstruction, FP64
updates, and solution slab writeback are included in solve timing. Common
input/operator preparation, shader compilation, and the subsequent independent
native verification are outside it. Hierarchy setup is reported separately.

GPU allocations at ready were 11.198 / 11.040 / 3.917 GiB for the three methods,
below the M1 Max's 24.960 GiB recommended working set. No numerical jobs
overlapped. The second basis and all cycle/solver parameters were held fixed.

[Machine-readable summary](results/near-null-svd24/summary.json),
[settings and input provenance](results/near-null-svd24/settings.json),
[new-MG Chroma checks](results/near-null-svd24/near.native.csv),
[filtered-MG checks](results/near-null-svd24/filtered.native.csv),
[BiCGStab checks](results/near-null-svd24/bicg.native.csv), and
[saved-field hashes](results/near-null-svd24/binary-artifacts.json) preserve
the complete comparison. The measured source is commit `f182d2a`, with
executable and per-source-file hashes recorded alongside the results.

## Construction

The fine operator is the retained-odd clover Schur complement

```
S = A_oo - H_oe A_ee^-1 H_eo.
```

For each extra test vector, the builder does the following:

1. Generate a deterministic dense complex Z4 vector from setup tags starting
   at 771337. These tags are independent of the physical solve sources.
2. Optionally apply damped clover-block Jacobi sweeps to the starting vector.
3. Form `b = S*y`, then approximately solve `S*z = b` from zero using GPU FP32
   FGMRES with local clover-block preconditioning. Form `v = y-z`. This is the
   homogeneous-solve error construction used in `mgproton`: a small `S*v`
   identifies a direction that an approximate inverse leaves behind.
4. Reapply `S` to measure the actual homogeneous residual, and normalize the
   whole candidate using a CPU FP64 norm. The iteration cap is allowed to
   return an inexact candidate; logs record actual reductions and iterations.
5. Gather candidates separately in each aggregate and each diagonal-gamma5
   chirality. Compute a complex thin SVD in CPU FP64 using Accelerate LAPACK
   `zgesvd`. Retain the leading **left** singular vectors as interpolation
   columns. With more candidates than retained vectors, this selects the
   locally most represented directions.
6. Check local FP64 orthogonality and numerical rank, then round interpolation
   weights to FP32 for application. Sparse Galerkin assembly remains FP64.

The entire oversampled candidate set is held in FP32. SVD workspace is local
to a worker, and the GPU setup Krylov workspace is released before hierarchy
application. This implementation uses one final SVD of all candidates, rather
than `mgproton`'s optional incremental, singular-value-weighted compression.
It implements full-vector candidate normalization before local compression.

The second interpolation retains the existing FP64 Jacobi construction, so
the comparison isolates the new fine-level basis. Both levels preserve
chirality. Every physical solve still uses FP64 defect corrections, full CPU
FP64 residual checks, and a separate original Chroma full-operator check.

## Options and checks

```
--basis filtered|near-null-svd
--null-vectors 24
--null-iters 128
--null-restart 12
--null-warm 32
--null-tol 0.001
```

The retained vector count is eight per chirality, matching the existing
low-memory driver. The positional filter-sweep argument applies to the
`filtered` mode; `--null-warm` controls the new mode's starting relaxation.
The null-solve target is a setup diagnostic, independent of physical solution
acceptance at `1e-12`.

Run the analytic complex SVD check without loading a gauge field or a GPU:

```sh
make -f Makefile.production build/low_memory_solve
VECLIB_MAXIMUM_THREADS=1 build/low_memory_solve --basis-self-test
```

It checks the complex singular spectrum `(3,2,1)`, optimal rank-two projection
with retained energy `13/14` in both chiralities, orthogonality, and rejection
of a rank-deficient candidate set. The separate real 12³×96 preflight uses two
physical RHS and the original
Chroma FP64 verifier; its evidence is in [near-null-svd12](results/near-null-svd12/).

For the 24³×256 light-mass comparison, use the same original Chroma bridge and
action metadata as the previous real-field study. The mapped low-memory driver
requires the aligned `CHRMG003` bridge layout; inserting its four-byte padding
does not change gauge or probe payload values. Exact input and executable
hashes are recorded in [near-null-svd24](results/near-null-svd24/).
Fresh native exports can write this layout directly with
`METALCLOVER_ALIGNED_BRIDGE=1`; see [LOW_MEMORY.md](LOW_MEMORY.md) for the
native exporter/verifier build requirements. Passing `-` for the exporter's
eigenvector input skips eigenvector extraction for these random-source tests.

```sh
VECLIB_MAXIMUM_THREADS=1 /private/tmp/mg-near-null-solve \
  /private/tmp/mg-null24 /private/tmp/mg-near-null24/near \
  mg 2 3,3,3,8 12 128 16 4 --basis near-null-svd \
  --null-vectors 24 --null-iters 128 --null-warm 32
```

Run the filtered control with the same positional arguments and
`--basis filtered`, and Metal BiCGStab with `bicg` in place of `mg`, using
fresh, distinct output prefixes; existing result CSVs are never overwritten.
Each process keeps one physical RHS resident and
reuses its hierarchy across both sources. Stage executables and binary fields
outside iCloud, and run numerical programs sequentially. The final saved
fields must pass `low_memory_verify` using the matching action XML and bridge.

After collecting the three timing CSVs and the native checks, regenerate the
validated summary with `python3 near_null_report.py results/near-null-svd24`.
The script rejects missing RHS, failed full residuals, incomplete runs, and
multiple hierarchy constructions.
