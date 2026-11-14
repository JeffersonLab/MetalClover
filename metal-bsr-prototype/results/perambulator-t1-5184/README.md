# Configuration 1000: standard perambulator at source/sink t=1

Repeat of the full 5,184-eigenvector calculation on the real 12³×96
configuration 1000, with four source spins (20,736 physical RHS). Both source
and retained sink times are one. The source is the supplied
`eigs_mod.n5184/1000/eigs.t_1.mod1000`, with its stored ordering and phases.
The anisotropic stout-clover fermion action is identical to the accepted
time-zero calculation; `action.xml` changes only `Contractions/t_sources`
from zero to one. Temporal boundary conditions remain antiperiodic.

The accepted point-sink prefix is
`/Users/admin/Downloads/MetalClover-results/cfg1000-t1-unsmeared`.
The standard perambulator filename will be
`/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_1.sdb1000`.
Its sixteen complex-FP64 spin-block matrices use source/sink time one in
the actual keys and metadata; `t0_1` in the filename labels source time one.

The production solve reuses one hierarchy for all 20,736 RHS in caller/device
batches of sixteen. Every solution must pass both the native qprop and the
independent full original Chroma CPU FP64 residual evaluations at relative
tolerance 1e-12 before retaining its time-one slice. Per-batch checkpoints
and complete column checksums support recovery and final verification.
The subsequent standard projection verifies input SHA256, direct QDP dot
products, complete Chroma value readback, full matrix gamma5 Hermiticity,
and independent DataLib reader compatibility, as at time zero.

`temporary-artifact-cleanup.json` records old regenerable benchmark outputs
removed from `/private/tmp` after matching their archived hashes. Previously
requested production outputs and all original gauge/eigenvector inputs are
retained. Code and small validation/provenance evidence are versioned; large
point-sink and perambulator files stay outside Git.


The preproduction validation passed 64 time-one MG RHS, a controlled
checkpoint/restart, and four matching BiCGStab RHS. The maximum full residuals
were 1.0360e-14 (MG) and 1.4814e-14 (BiCGStab), with point-slice relative
solution difference 6.8601e-14. The 16-vector projected pilot had gamma5
Hermiticity relative Frobenius discrepancy 6.6575e-16; all 4,096 complex
values matched bitwise through the independent DataLib analysis reader.
Wrong eigenvector time metadata was rejected before output creation. The
new time-default-zero executable verified an existing complete time-zero
checkpoint, and all sixteen value hashes of the time-zero projection
regression were unchanged. The interrupted pilot intentionally builds two
hierarchies across its two sessions; production starts fresh and uses one.
