# Native Chroma adapter validation

The final source builds on Chroma `rge/devel-datalib-reorg` base e43d28927
with the installable MetalClover C ABI. `validation.json` records configuration,
input/source hashes, eight whole-lattice operator probes and all 20 accepted
original-Chroma FP64 scalar, reconstruction, qprop-array and actual
`SB::ChimeraSolver`/`SB::doInversion` checks. Each MG object retains its one
hierarchy across two RHS. The array distillation check uses two deterministic
complex colorvectors, source spin 3, t=80, and `max_rhs=1`.

The test uses the real 12³×96 field and its own action parameters (mass −0.0743).
It is separate from the standalone 32³×256 mass −0.0860 performance comparison.
Full Chroma 32-cubed memory fit and native integration performance have not
been established. The native harness uses additional in-memory FP64 fields.

The enabled full build and disabled factory-aggregate compile passed. Logs
include pre-existing macro, override and duplicate-library warnings; the full
build succeeded. Large full build logs remain local and their paths/hashes
are recorded in `build-log-manifest.json`.
