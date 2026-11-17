# Configuration 1000: standard perambulator at source/sink t=6

**Completed and verified, 2026-11-16 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,202.497 s (8 h 23 min 22 s) |
| Native qprop calls, including native full FP64 checks | 28,025.631 s |
| Additional independent full FP64 checks | 1,967.618 s |
| Point extraction/write/readback | 31.636 s |
| Maximum full relative residual, either check | 2.434164278807051e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.954 GiB |
| Nominal dispatched GPU operator work | 8,388,284,141,297,664 FLOP |
| Nominal GPU operator work / total solve wall | 277.73 GFLOP/s |
| Standard projection plus complete native validation | 274.779 s |
| Complex-double GEMM time | 35.454 s |
| Maximum sampled projection process RSS | 11.368 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.09839082579092e-14 |
| Maximum absolute gamma5 discrepancy | 1.38751204141256e-15 |
| Maximum independent QDP contraction error, norm-scaled | 4.76180028260775e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,504.906 s (8 h 28 min 25 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_6.sdb1000
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
6 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`78f6f1a613c67a3f55cdc9e271de644bb94e8f22b1359dae933afe5e92b69b37`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1000-t6-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`fb5e42afd07c6fa3966399d0ff3eb3f6e25dc5bc6fc766e00e2eefb1982cd43c`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_6.mod1000`,
215,332,759 bytes, SHA256
`8379a7caa1e2c5b2b44222e8c99dea7e02104777f488ab6ac67a958b28ad19c8`. Their header reports 12³,
decay_dir=3, t_slice=6 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 6 in `action.xml`. The mass −0.0743,
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
and complete logs/residual evidence remain. See [the retirement record](../perambulator-t6-5184/intermediate-retirement.json).
