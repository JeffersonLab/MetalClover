# MetalClover

Standalone Apple GPU experiments for anisotropic clover-Wilson inversion.

## Provenance

This project grew out of studying Eloy Romero's `mgproton` multigrid prototype
in Chroma and his [SuperbBLAS](https://github.com/eromero-vlc/superbblas) work.
We credit that work for guiding the initial sparse color/spin representation
and the investigation of even/odd preconditioning, chirality-split aggregation,
Galerkin coarse operators, and Krylov solvers. MetalClover implements these
ideas with dedicated CPU and Metal kernels. The optional near-null/SVD basis
builder follows the homogeneous-solve and local-compression approach in
`mgproton`; the hierarchy application and memory management remain separate
implementations.

The fermion-action conventions and independent numerical references come from
[Chroma](https://github.com/JeffersonLab/chroma) and QDP++. The anisotropic
stout-smeared clover-Wilson operator follows Chroma's conventions, and the native
adapter and verification tools use Chroma/QDP++ to check the stout gauge state,
operator application, and full FP64 residuals. These checks establish agreement
with the Chroma operator; they do not establish identical hierarchy construction
or convergence behavior to `mgproton`. We acknowledge Eloy Romero and the
Chroma/QDP++ developers for these foundations.

## Code and results

The source and reproducible numerical results are in [metal-bsr-prototype](metal-bsr-prototype/).

- [Near-null/SVD MG on 32³×256: 4.61× warm solve speedup and original Chroma FP64 checks](metal-bsr-prototype/NEAR_NULL_SVD32.md)
- [Real 32³×256 low-memory LDL clover, tuned MG timings, and native FP64 checks](metal-bsr-prototype/LOW_MEMORY.md)
- [Real 24³×256 MG tuning and original Chroma FP64 validation](metal-bsr-prototype/MG_TUNING_REAL.md)
- [Homogeneous near-null solves and SVD-compressed MG basis](metal-bsr-prototype/NEAR_NULL_SVD.md)
- [Native Chroma adapter, action conventions, and initial real-field timings](metal-bsr-prototype/PERFORMANCE_REAL_FIELDS.md)
- [Standalone prototype overview](metal-bsr-prototype/README.md)

Large gauge fields, generated binaries, and saved solution fields are excluded from Git;
text results, input/output hashes, source, and reproduction instructions are versioned.

## Library and Chroma integration

Build the installable Metal library with Apple Clang:

```sh
make install PREFIX=/path/to/metal-prefix
```

All C++ library types and routines are in `MetalCloverEnv`. The public header
`metalclover/metalclover.h` provides a C ABI usable from GCC-built Chroma,
separate full and odd-Schur solve interfaces, arrays of source/solution pointers,
and structured residual, iteration, timing and nominal operator-work results.
One context owns and reuses one hierarchy. Array solves run batched GPU tiles
(default cap 8), sharing gauge, clover and hierarchy data. The caller owns
the source and solution arrays.
See [batched RHS support](metal-bsr-prototype/BATCH_RHS.md) and
[interfaces and FLOP accounting](metal-bsr-prototype/OPERATOR_WORK.md).

The Chroma integration on `rge/devel-datalib-reorg` enables
`MetalBiCGClover` and `MetalMGClover` with `--with-metal=/path/to/metal-prefix`.
It implements the vector SystemSolver interface used by the SuperbBLAS
propagator/distillation measurement, including native even/odd reduction and
reconstruction. Each refinement uses the supplied original Chroma FP64
operator. See Chroma's `docs/metal-clover.md` for XML and the real-gauge harness.

## Formatting and continuous integration

Run `make format` from the repo root or `metal-bsr-prototype/` to format the
`.hpp`, `.cc`, `.mm`, and `.metal` sources using the checked-in
[`.clang-format`](.clang-format). This requires `clang-format` on `PATH`;
override the executable with `make format CLANG_FORMAT=clang-format-23` if needed.
Use `make format-check` to check formatting without modifying files. CI pins
clang-format to **23.1.5**, the version used to format the current sources.

[GitHub Actions](.github/workflows/ci.yml) runs on pushes, pull requests, and
manual dispatches. It checks formatting on Ubuntu and builds all 19 standalone
Objective-C++ programs on an Apple Silicon `macos-15` runner with Xcode 16.4.
It also compiles all ten Metal shader files offline using `make -C
metal-bsr-prototype metal-check`; this does not require GPU execution.

CI covers the standalone build. Native Chroma/QDP++ adapters require a separate
installation, and the real-field residual checks and performance studies require
local gauge data and an Apple GPU. Those numerical runs remain separate from
the hosted formatting and build checks.
