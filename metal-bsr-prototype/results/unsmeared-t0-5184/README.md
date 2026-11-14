# Configuration 1000: all 5,184 eigenvector sources, unsmeared t=0 sinks

The requested production collection contains all four source spins for each
supplied t=0 eigenvector: 20,736 physical RHS. Source and sink time are both
zero; sinks retain the full unsmeared site/color/spin slice. See
[the writer guide](../../UNSMEARED_T0.md) for axes, native Chroma checks and
restart instructions. Large output files stay outside Git.

**Completed and verified on 2026-11-14.** All 20,736 columns passed both full
four-dimensional original Chroma CPU FP64 residual checks at 1e-12. The native
completed-output verifier then reread every column checksum successfully;
the controller validated the complete NPY header, size, source/spin ordering,
norms and residual journal and computed the whole-file SHA256.

| Measurement | Result |
| --- | ---: |
| Physical RHS | 20,736 |
| Caller/device RHS cap | 16 / 16 |
| Array calls | 1,296 |
| Hierarchies built / sessions | 1 / 1 |
| End-to-end solve-session wall time | 30,140.426 s = 8 h 22 min 20 s |
| Wall seconds per physical RHS | 1.45353 |
| Native qprop calls, including native FP64 checks | 27,968.787 s |
| Additional independent full CPU FP64 checks | 1,963.553 s |
| Output extraction/write/readback | 28.886 s |
| Solver/hierarchy construction | 21.041 s |
| Iterations per physical RHS | 73–82 |
| Maximum full relative residual, either check | 2.4063810813519627e-14 |
| Peak process RSS | 5.955 GiB |
| Nominal dispatched GPU operator work | 8,363,673,104,130,048 FLOP |
| Nominal operator rate over total session wall | 277.49 GFLOP/s |

The machine was an Apple M1 Max with 24 GPU cores and 32 GiB system memory,
using one rank. The wall time includes source loading, setup, source/solution
copies, all residual checks, output and checkpointing; it excludes the final
completed-file checksum reread and SHA256. The operator-work numerator counts
dispatched GPU operators, including inactive/padded lanes. It excludes CPU
FP64, Krylov-vector arithmetic and hierarchy setup work. This rate uses the
same nominal-work model as the inverter's returned statistics.

Each source spin has 5,184 accepted columns. Maximum independent full FP64
residuals for source spins 0, 1, 2 and 3 were respectively 2.4064e-14,
2.3882e-14, 2.3848e-14 and 2.3541e-14. Source norms ranged from
0.9999999986416703 to 1.0000000014039585. These checks apply to the supplied
FP32 eigenvectors promoted to FP64, without changing their phases or order.

The production prefix is
`/Users/admin/Downloads/MetalClover-results/cfg1000-t0-unsmeared`.
Its `.run.json` reports live progress; `.complete.json` appears only after
completion, native column checksum rereading, journal validation and SHA256.
The source is the supplied descending-order complex-FP32 MOD basis, promoted
to FP64 without sorting or phase rotation.

`action.xml` records the original anisotropic stout-clover parameters,
including mass −0.0743 and antiperiodic temporal boundary conditions.
One MetalMGClover solver and hierarchy span an uninterrupted production run.
Caller and GPU caps are sixteen, with source spin as the outer loop.
Every accepted column must pass two full four-dimensional original Chroma
CPU FP64 residual checks at relative tolerance 1e-12.

The retained `.npy` file is 6,879,707,264 bytes, complex128, with shape
`(5184,4,12,12,12,3,4)`. Its SHA256 is
`3b017ee3fd6989514f3ac23f79859082442da18b3fb201c278d14fd85c8ffc51`.
This is the unsmeared time-zero sink block with eigenvector/spin sources;
the full four-dimensional solutions were checked before discarding other
sink times.

The `cfg1000-t0-unsmeared.*` text sidecars contain the complete journal,
metadata, native log, resource samples, completion report and final verification
log. `production-validation.json` checks all four source-spin collections,
all 1,296 width-sixteen reused-hierarchy calls, the single setup, and unchanged
frozen source/executable/library hashes. The production controller is commit
`60e62ed53915ce464c107a9ab82e2827b620a355`; its Chroma writer is commit
`6aeb610a133b83ec910e0530471c34494117fde1`.

`validation/pilot-summary.json` records actual highest/lowest-order basis
sources, four-spin output, native and independent FP64 acceptance, point-sink
MG/BiCGStab agreement, per-column payload hashes, and recovery from deliberately
uncommitted test data after a checkpoint. The final-width pilot separately
checks the complete controller and corrected nested XML metadata. Frozen
source, executable and library SHA256 hashes are recorded alongside it.
Numerical solver/library source was unchanged for this writer addition.
