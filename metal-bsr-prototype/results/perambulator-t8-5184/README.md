# Configuration 1000: standard perambulator at source/sink t=8

**Completed and verified, 2026-11-17 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,178.227 s (8 h 22 min 58 s) |
| Native qprop calls, including native full FP64 checks | 28,008.504 s |
| Additional independent full FP64 checks | 1,964.601 s |
| Point extraction/write/readback | 31.494 s |
| Maximum full relative residual, either check | 2.447499776549151e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.960 GiB |
| Nominal dispatched GPU operator work | 8,396,159,228,805,120 FLOP |
| Nominal GPU operator work / total solve wall | 278.22 GFLOP/s |
| Standard projection plus complete native validation | 273.023 s |
| Complex-double GEMM time | 35.519 s |
| Maximum sampled projection process RSS | 11.368 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.1185654441594e-14 |
| Maximum absolute gamma5 discrepancy | 1.45973200563691e-15 |
| Maximum independent QDP contraction error, norm-scaled | 4.77226705234107e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,479.021 s (8 h 27 min 59 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_8.sdb1000
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
8 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`7862a06abe5b68b3b62616c19f6975e1109280b165bf05e5915759f318fed9b6`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1000-t8-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`ca5d5b85364f6e9cec991aa7bbb60da661e13aafbddbbcd5216f97d1bfb5be8e`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_8.mod1000`,
215,332,808 bytes, SHA256
`7b794f08efa38d2fca78f168db7ccb5bbf18e19b3e4dbb0715aca371d86ca0e7`. Their header reports 12³,
decay_dir=3, t_slice=8 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 8 in `action.xml`. The mass −0.0743,
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
and complete logs/residual evidence remain. See [the retirement record](../perambulator-t8-5184/intermediate-retirement.json).
