# Configuration 1000: standard perambulator at source/sink t=9

**Completed and verified, 2026-11-17 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 31,557.584 s (8 h 45 min 58 s) |
| Native qprop calls, including native full FP64 checks | 29,260.316 s |
| Additional independent full FP64 checks | 2,080.640 s |
| Point extraction/write/readback | 25.993 s |
| Maximum full relative residual, either check | 2.467216318862918e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.959 GiB |
| Nominal dispatched GPU operator work | 8,403,223,803,322,368 FLOP |
| Nominal GPU operator work / total solve wall | 266.28 GFLOP/s |
| Standard projection plus complete native validation | 279.950 s |
| Complex-double GEMM time | 36.994 s |
| Maximum sampled projection process RSS | 11.337 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.13302490016386e-14 |
| Maximum absolute gamma5 discrepancy | 1.55695563358594e-15 |
| Maximum independent QDP contraction error, norm-scaled | 2.74147461829395e-15 |
| Complete local pipeline, including whole-file hashes and reader | 31,866.067 s (8 h 51 min 6 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_9.sdb1000
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
9 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`18f98412e66a21087dd681e96843e0268e7719aa7be34b7986b77b9b99880174`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1000-t9-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`d61abaf5695f5c16eb40f26d2ad2332e1c729fc02e30ae0ef10abcc6c9a17d3d`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_9.mod1000`,
215,332,765 bytes, SHA256
`5d2e7d8b02a510cb9efb83371d3a33101a63a3c1962083c600f1b9d957df6e95`. Their header reports 12³,
decay_dir=3, t_slice=9 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 9 in `action.xml`. The mass −0.0743,
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
and complete logs/residual evidence remain. See [the retirement record](../perambulator-t9-5184/intermediate-retirement.json).
