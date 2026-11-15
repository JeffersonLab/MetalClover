# Configuration 1000: standard perambulator at source/sink t=3

**Completed and verified, 2026-11-15 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,172.532 s (8 h 22 min 53 s) |
| Native qprop calls, including native full FP64 checks | 28,001.372 s |
| Additional independent full FP64 checks | 1,964.159 s |
| Point extraction/write/readback | 31.342 s |
| Maximum full relative residual, either check | 2.464448516586995e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.954 GiB |
| Nominal dispatched GPU operator work | 8,385,719,574,061,056 FLOP |
| Nominal GPU operator work / total solve wall | 277.93 GFLOP/s |
| Standard projection plus complete native validation | 285.502 s |
| Complex-double GEMM time | 35.642 s |
| Maximum sampled projection process RSS | 11.368 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.08170893841624e-14 |
| Maximum absolute gamma5 discrepancy | 1.47832923494512e-15 |
| Maximum independent QDP contraction error, norm-scaled | 3.46972980107725e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,482.898 s (8 h 28 min 3 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_3.sdb1000
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
3 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`346a37dbea6ed810e8a209e42fc482c895366024dc12e00d199466731e88ebd8`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1000-t3-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`47ea5a8dbe16f7a72ec9f93dd698eaaaf684c26722210ec6d1c391402a166d5d`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_3.mod1000`,
215,332,763 bytes, SHA256
`b2b5b12df39c2952c3473e59ae6fb9e4decfd1b8ea6ef853b00506dbbdd1a79e`. Their header reports 12³,
decay_dir=3, t_slice=3 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 3 in `action.xml`. The mass −0.0743,
anisotropies, clover coefficients, spatial stout action and antiperiodic
temporal boundary conditions are identical to the earlier accepted runs.
No numerical source, executable, library or solver setting changed.

`production-validation.json` checks every source-spin collection, all native
array calls, one hierarchy/session, both residual bounds, complete file checks
and frozen source/executable/library hashes. `point.*` and `peram.*` contain
the full journals, metadata, logs, checksums and resource samples. The queue
and pipeline scripts run locally with sleep prevention, independently of the
user's iPhone connection. Large output data remain outside Git.

Storage: with explicit user approval, the earlier time-zero and time-one
unsmeared `.npy` intermediates were removed after checking their hashes and
rechecking the retained standard perambulator hashes. Every standard file,
original gauge/eigenvector input, log and residual record is retained; see
[the cleanup record](../perambulator-t3-5184/intermediate-cleanup.json).


Storage update for the time-four repeat: the accepted time-two and time-three
unsmeared `.npy` intermediates were retired after all checks and evidence capture,
rechecking both point and retained standard SHA256 hashes, and confirming the
intermediates were not open. Every standard perambulator, original input, log
and accepted residual record is retained. See
[the time-four cleanup record](../perambulator-t4-5184/intermediate-cleanup.json).
