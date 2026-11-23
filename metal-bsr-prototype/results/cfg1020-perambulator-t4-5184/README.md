# Configuration 1020: standard perambulator at source/sink t=4

**Completed and verified, 2026-11-19 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 29,953.624 s (8 h 19 min 14 s) |
| Native qprop calls, including native full FP64 checks | 27,784.721 s |
| Additional independent full FP64 checks | 1,962.937 s |
| Point extraction/write/readback | 31.626 s |
| Maximum full relative residual, either check | 2.487838461012357e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.953 GiB |
| Nominal dispatched GPU operator work | 8,347,196,466,511,872 FLOP |
| Nominal GPU operator work / total solve wall | 278.67 GFLOP/s |
| Standard projection plus complete native validation | 275.199 s |
| Complex-double GEMM time | 35.581 s |
| Maximum sampled projection process RSS | 11.368 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.09986128214969e-14 |
| Maximum absolute gamma5 discrepancy | 1.52355268089469e-15 |
| Maximum independent QDP contraction error, norm-scaled | 5.47264608052771e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,254.593 s (8 h 24 min 15 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_4.sdb1020
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
4 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`c17085161ddcb12caf5771de6c32817a025a9a7b7ec326f1d92644b57572fb5f`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1020-t4-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`0dbfdca4ddfbbab2da4ec5b22fe9e8419f7a36b7fb41a5c266c24c1ca1b91d8c`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1020/eigs.t_4.mod1020`,
215,332,809 bytes, SHA256
`661dded1b4479a2d8c1f0983ccbe34be4913e34acf05395aaeabec238307aa45`. Their header reports 12³,
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
and complete logs/residual evidence remain. See [the retirement record](../cfg1020-perambulator-t4-5184/intermediate-retirement.json).
