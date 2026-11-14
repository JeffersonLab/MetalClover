# Restored Eloy spin-major caller with Metal batching

The production perambulator measurement and header, and the helper header,
match their pre-overhaul Chroma `b3a380ea6` contents byte-for-byte. The helper
implementation restores the original source-spin outer loop with eigenvector
batches inside; useful physical-RHS logging is retained. Localized bounded
temporary pools, the `Ns` spin dimension and optional PRIMME guard remain.
The Metal array backend and numerical dylib are unchanged.

On the real 12³×96 LIME field at mass −0.0743, source time zero, all 64 native
checks passed. For both methods, one vector at `max_rhs=1` produced widths
1,1,1,1. Three vectors at `max_rhs=2` produced 2,1,2,1,2,1,2,1. Full original
Chroma CPU FP64 relative residuals were required to be at most 1e-12, and
selected source-time contractions were compared against direct CPU inner
products. The 3+2 partial-spin path, scalar/odd/full/reconstructed paths,
nonzero/exact guesses and zero RHS checks remain included.

The spin-major harness smoke run compares caller caps 1/8/16 on one vector
(four inversions) and then two vectors (eight). Actual caller/GPU widths,
four-spin accounting, full-field scalar/BiCGStab/MG agreement and one reused
MG hierarchy were verified. These small timings are not a fresh 96/192-RHS
performance measurement. The full inline perambulator file-I/O job and the
full 5,184-vector job were not run.

`summary.json` contains acceptance and observed batches; `sources.json` hashes
the exact canonical/tested source copies. `inputs.json` hashes the gauge and
action XML. `artifacts.json` records binaries, the unchanged Metal dylib and
local scalar reference file; those binaries and reference fields are excluded
from Git. Build logs, reproduction controllers and normalized log hashes are
included. Copied text logs use LF and trim trailing whitespace.

The matching Chroma restoration is
[`d8765a055`](https://github.com/JeffersonLab/chroma/commit/d8765a05538285412b7983888ab336546ea0009f)
on `rge/devel-datalib-reorg`. The MetalClover documentation/evidence update is
the commit containing this directory. `integration.json` records the pairing.
The restoration script records the edit against the preceding coalescing
version; it is not an idempotent patch for an already restored checkout.
