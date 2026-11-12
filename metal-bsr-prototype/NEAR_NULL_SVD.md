# Homogeneous near-null solves and chiral SVD compression

`low_memory_solve` now has an optional `--basis near-null-svd` builder inspired
by Eloy Romero's `getMGProlongator` in `mgproton.cc`. The original filtered
builder remains the default and is selected explicitly by `--basis filtered`.
The fine action, Schur operator, local LDL factors, physical RHS, outer
refinement, and original Chroma verifier are shared by both modes.

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

It checks a known complex singular spectrum, optimal rank-two projection in
both chiralities, orthogonality, and rejection of a rank-deficient candidate
set. The separate real 12³×96 preflight uses two physical RHS and the original
Chroma FP64 verifier; its evidence is in [near-null-svd12](results/near-null-svd12/).

For the 24³×256 light-mass comparison, use the same original Chroma bridge and
action metadata as the previous real-field study. The mapped low-memory driver
requires the aligned `CHRMG003` bridge layout; inserting its four-byte padding
does not change gauge or probe payload values. Exact input and executable
hashes are recorded in [near-null-svd24](results/near-null-svd24/).

```sh
VECLIB_MAXIMUM_THREADS=1 /private/tmp/mg-near-null-solve \
  /private/tmp/mg-null24 /private/tmp/mg-near-null24/near \
  mg 2 3,3,3,8 12 128 16 4 --basis near-null-svd \
  --null-vectors 24 --null-iters 128 --null-warm 32
```

Run the filtered control with the same positional arguments and
`--basis filtered`, and Metal BiCGStab with `bicg` in place of `mg`, using
distinct output prefixes. Each process keeps one physical RHS resident and
reuses its hierarchy across both sources. Stage executables and binary fields
outside iCloud, and run numerical programs sequentially. The final saved
fields must pass `low_memory_verify` using the matching action XML and bridge.
