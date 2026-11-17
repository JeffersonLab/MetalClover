# Configuration 1000: standard perambulator at source/sink t=4

**Completed and verified, 2026-11-16 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,316.207 s (8 h 25 min 16 s) |
| Native qprop calls, including native full FP64 checks | 28,132.907 s |
| Additional independent full FP64 checks | 1,974.947 s |
| Point extraction/write/readback | 33.094 s |
| Maximum full relative residual, either check | 2.405119835565921e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.954 GiB |
| Nominal dispatched GPU operator work | 8,387,569,482,596,352 FLOP |
| Nominal GPU operator work / total solve wall | 276.67 GFLOP/s |
| Standard projection plus complete native validation | 271.954 s |
| Complex-double GEMM time | 35.380 s |
| Maximum sampled projection process RSS | 11.375 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.07964511896396e-14 |
| Maximum absolute gamma5 discrepancy | 1.59575750130691e-15 |
| Maximum independent QDP contraction error, norm-scaled | 4.76936098257876e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,617.778 s (8 h 30 min 18 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_4.sdb1000
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
4 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`9f9151a3c6b52c22016a43948045089a4c5f4d9770e03c9bac933b95757abf23`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1000-t4-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`6c674c8f02f30755309cdb4adbb8d310aef7c276e5a1e76d912c2f183e414a57`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_4.mod1000`,
215,332,828 bytes, SHA256
`5605ffb6568cc8fd3800c66cf09bbe21be740abdf254612f73624261372a6114`. Their header reports 12³,
decay_dir=3, t_slice=4 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 4 in `action.xml`. The mass −0.0743,
anisotropies, clover coefficients, spatial stout action and antiperiodic
temporal boundary conditions are identical to the earlier accepted runs.
No numerical source, executable, library or solver setting changed.

`production-validation.json` checks every source-spin collection, all native
array calls, one hierarchy/session, both residual bounds, complete file checks
and frozen source/executable/library hashes. `point.*` and `peram.*` contain
the full journals, metadata, logs, checksums and resource samples. The queue
and pipeline scripts run locally with sleep prevention, independently of the
user's iPhone connection. Large output data remain outside Git.

The accepted point intermediate was retired after complete verification and
point/standard SHA256 rechecks. Its standard perambulator, all original inputs
and complete logs/residual evidence remain. See [the retirement record](../perambulator-t6-5184/dependency-retirement.json).
