# Configuration 1020: standard perambulator at source/sink t=5

**Completed and verified, 2026-11-20 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,057.439 s (8 h 20 min 57 s) |
| Native qprop calls, including native full FP64 checks | 27,887.567 s |
| Additional independent full FP64 checks | 1,964.988 s |
| Point extraction/write/readback | 33.644 s |
| Maximum full relative residual, either check | 2.335408712697519e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.954 GiB |
| Nominal dispatched GPU operator work | 8,351,638,967,672,832 FLOP |
| Nominal GPU operator work / total solve wall | 277.86 GFLOP/s |
| Standard projection plus complete native validation | 279.424 s |
| Complex-double GEMM time | 35.556 s |
| Maximum sampled projection process RSS | 11.376 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.09406558553181e-14 |
| Maximum absolute gamma5 discrepancy | 1.45448701651045e-15 |
| Maximum independent QDP contraction error, norm-scaled | 3.31901541622502e-15 |
| Complete local pipeline, including whole-file hashes and reader | 30,365.372 s (8 h 26 min 5 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_5.sdb1020
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
5 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`9ff68602a9447adf6eb2f38d039102b17b8eb83cb85797a376873d1644288483`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1020-t5-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`1843bf0a35e9b94f839bfcb8dc2b859a889efe73ba416eda9a92875e4261c5b3`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1020/eigs.t_5.mod1020`,
215,332,789 bytes, SHA256
`9ba7649b3585670b145d68e7345a76f6598806243a30738f43352f13e4241197`. Their header reports 12³,
decay_dir=3, t_slice=5 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 5 in `action.xml`. The mass −0.0743,
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
and complete logs/residual evidence remain. See [the retirement record](../cfg1020-perambulator-t5-5184/intermediate-retirement.json).
