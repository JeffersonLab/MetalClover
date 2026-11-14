# Configuration 1000: all 5,184 eigenvector sources, unsmeared t=0 sinks

The requested production collection contains all four source spins for each
supplied t=0 eigenvector: 20,736 physical RHS. Source and sink time are both
zero; sinks retain the full unsmeared site/color/spin slice. See
[the writer guide](../../UNSMEARED_T0.md) for axes, native Chroma checks and
restart instructions. Large output files stay outside Git.

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

`validation/pilot-summary.json` records actual highest/lowest-order basis
sources, four-spin output, native and independent FP64 acceptance, point-sink
MG/BiCGStab agreement, per-column payload hashes, and recovery from deliberately
uncommitted test data after a checkpoint. The final-width pilot separately
checks the complete controller and corrected nested XML metadata. Frozen
source, executable and library SHA256 hashes are recorded alongside it.
Numerical solver/library source was unchanged for this writer addition.
