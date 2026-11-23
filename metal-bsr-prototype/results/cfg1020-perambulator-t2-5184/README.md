# Configuration 1020: standard perambulator at source/sink t=2

**Completed and verified, 2026-11-18 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,074.512 s (8 h 21 min 15 s) |
| Native qprop calls, including native full FP64 checks | 27,910.041 s |
| Additional independent full FP64 checks | 1,963.716 s |
| Point extraction/write/readback | 29.566 s |
| Maximum full relative residual, either check | 2.571971200775573e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.954 GiB |
| Nominal dispatched GPU operator work | 8,359,163,858,681,856 FLOP |
| Nominal GPU operator work / total solve wall | 277.95 GFLOP/s |
| Standard projection plus complete native validation | 275.073 s |
| Complex-double GEMM time | 35.534 s |
| Maximum sampled projection process RSS | 11.376 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.11850094889157e-14 |
| Maximum absolute gamma5 discrepancy | 1.46485560070912e-15 |
| Maximum independent QDP contraction error, norm-scaled | 4.06483780712527e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,375.766 s (8 h 26 min 16 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_2.sdb1020
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
2 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`e5642ec7f3c25d24e4638f50ad4ca92ba28dc8dff319c9b695b80994a76b4047`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1020-t2-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`add38d294ffbe3751a41495c1cf94c75a1d4c64526284fcf6a871471bdac61e1`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1020/eigs.t_2.mod1020`,
215,332,772 bytes, SHA256
`742e34d2d9f049a69241e2c5b68f4f8c6129d2c7ad68e7d4af4c1c556abbcac2`. Their header reports 12³,
decay_dir=3, t_slice=2 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 2 in `action.xml`. The mass −0.0743,
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
and complete logs/residual evidence remain. See [the retirement record](../cfg1020-perambulator-t2-5184/intermediate-retirement.json).
