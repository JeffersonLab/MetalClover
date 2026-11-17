# MetalClover documentation

- [API guide](API.md): public C ABI, data layout, options, full and odd solves,
  multi-RHS tiles, callbacks, errors, and a C++ usage example.
- [Build and checks](BUILD.md): library installation, client linking, standalone
  programs, shader checks, numerical tests, and optional Chroma integration.
- [GPU implementation notes (LaTeX)](GPU_IMPLEMENTATION_NOTES.tex): Schur
  multigrid mathematics, Galerkin operators, transfers, near-null/SVD setup,
  Apple Metal implementation choices, and mgproton provenance.

The public contract is declared in
[`metalclover.h`](../metal-bsr-prototype/include/metalclover/metalclover.h).
The internal C++ headers and benchmark drivers are experimental implementation
interfaces; use the installed C ABI for application integration.

Detailed historical studies remain under `metal-bsr-prototype/`, including
[multi-RHS measurements](../metal-bsr-prototype/BATCH_RHS.md),
[distillation source conventions](../metal-bsr-prototype/DISTILLATION_RHS.md),
[operator-work accounting](../metal-bsr-prototype/OPERATOR_WORK.md), and
[large-volume SVD results](../metal-bsr-prototype/NEAR_NULL_SVD32.md).

Paths named in the LaTeX source refer to files under `metal-bsr-prototype/`
unless explicitly described as root files. To compile the notes from the root:

```sh
mkdir -p /tmp/metalclover-notes
pdflatex -halt-on-error -output-directory=/tmp/metalclover-notes docs/GPU_IMPLEMENTATION_NOTES.tex
pdflatex -halt-on-error -output-directory=/tmp/metalclover-notes docs/GPU_IMPLEMENTATION_NOTES.tex
```
