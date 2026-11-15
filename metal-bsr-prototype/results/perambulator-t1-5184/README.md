# Configuration 1000: standard perambulator at source/sink t=1


**Completed and verified, 2026-11-15 UTC.** All 20,736 RHS passed both
full original Chroma CPU FP64 residual checks at 1e-12. One hierarchy served
all 1,296 width-sixteen array calls in one uninterrupted session. The complete
point payload passed native column checksum rereading and SHA256 verification.
The standard perambulator then passed all sixteen full record readbacks,
direct QDP contractions, full gamma5 Hermiticity and DataLib key checks.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,095.368 s (8 h 21 min 35 s) |
| Native qprop calls, including native full FP64 checks | 27,925.926 s |
| Additional independent full FP64 checks | 1,964.694 s |
| Point extraction/write/readback | 33.609 s |
| Maximum full relative residual, either check | 2.3818484661997255e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.954 GiB |
| Nominal dispatched GPU operator work | 8,368,542,173,995,008 FLOP |
| Nominal GPU operator work / total solve wall | 278.07 GFLOP/s |
| Standard projection plus complete native validation | 273.205 s |
| Complex-double GEMM time | 35.303 s |
| Maximum sampled projection process RSS | 11.372 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.08109966010116e-14 |
| Maximum absolute gamma5 discrepancy | 1.77129036300412e-15 |
| Maximum independent QDP contraction error, norm-scaled | 5.73857748774942e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,397.546 s (8 h 26 min 38 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator FLOP numerator includes dispatched inactive/padded lanes and excludes
CPU FP64 and Krylov-vector arithmetic; its denominator includes all solve-session
work. Final whole-file verification/hashing and standard projection are separate
from that solve rate. Process RSS is sampled every five seconds.

Point output SHA256:
`a45f5451781bc8b2ad58f22cc16c026408dfc502213f843cfaabf55539087a27`.
Standard perambulator SHA256:
`820f87e336648a40cd579043ba87efa6e03dd5fa4e935179e16f83778dd88cc5`.
The point file is 6,879,707,264 bytes; the standard file is 6,907,158,528 bytes.
`point.*` and `peram.*` record full journals, metadata, logs, verification and
resource samples. `production-validation.json` checks all four spin collections,
all native calls, unchanged frozen source/binary/library hashes and final acceptance.
The original input eigenvectors remain FP32; the accuracy checks refer to those
supplied vectors promoted to FP64, with order and phases preserved.

The successful time-one pipeline triggered the queued time-two run automatically.
Time-two uses the identical solver code and action, changing only the source/sink
slice and its supplied eigenvectors. Its independent hierarchy is built once
for that new full collection.

Repeat of the full 5,184-eigenvector calculation on the real 12³×96
configuration 1000, with four source spins (20,736 physical RHS). Both source
and retained sink times are one. The source is the supplied
`eigs_mod.n5184/1000/eigs.t_1.mod1000`, with its stored ordering and phases.
The anisotropic stout-clover fermion action is identical to the accepted
time-zero calculation; `action.xml` changes only `Contractions/t_sources`
from zero to one. Temporal boundary conditions remain antiperiodic.

The accepted point-sink prefix is
`/Users/admin/Downloads/MetalClover-results/cfg1000-t1-unsmeared`.
The verified standard perambulator is saved at
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


Storage update for the queued time-three run: with explicit user approval, the
intermediate time-zero and time-one unsmeared `.npy` payloads were removed after
checking their hashes and the retained standard perambulator hashes. All standard
`.sdb1000` files and the original inputs, logs and checksums remain. See
[the cleanup record](../perambulator-t3-5184/intermediate-cleanup.json).
