# Configuration 1020: standard perambulator at source/sink t=6

**Completed and verified, 2026-11-23 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 30,090.162 s (8 h 21 min 30 s) |
| Native qprop calls, including native full FP64 checks | 27,914.550 s |
| Additional independent full FP64 checks | 1,966.884 s |
| Point extraction/write/readback | 32.505 s |
| Maximum full relative residual, either check | 2.586792573545587e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.955 GiB |
| Nominal dispatched GPU operator work | 8,340,964,774,772,736 FLOP |
| Nominal GPU operator work / total solve wall | 277.20 GFLOP/s |
| Standard projection plus complete native validation | 277.561 s |
| Complex-double GEMM time | 37.330 s |
| Maximum sampled projection process RSS | 11.337 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.08440794058118e-14 |
| Maximum absolute gamma5 discrepancy | 1.35394028223089e-15 |
| Maximum independent QDP contraction error, norm-scaled | 4.46724965601697e-15 |
| Complete local pipeline, including whole-file hashes and reader | 337,495.374 s (93 h 44 min 55 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_6.sdb1020
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
6 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`a86a574a22e6195d369003d63d4658ebce971ef8373987b7af5af146d8ac0cdf`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1020-t6-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`6fa8aa2af61feb5ffe8425150e4e4f942b9cd4f8e041e50b0e8d0b8213201163`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1020/eigs.t_6.mod1020`,
215,332,791 bytes, SHA256
`453befe3ada05360190239c3a8736659ff0fee2a3d90cee7851d80b032fd4104`. Their header reports 12³,
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
and complete logs/residual evidence remain. See [the retirement record](../cfg1020-perambulator-t6-5184/intermediate-retirement.json).


## Postprocessing recovery

The native solve finished successfully on November 20. Missing temporary runtime
launchers prevented its subsequent checksum verifier from starting, followed by
the user-requested November 21 disk-space pause. The user resumed on November 23.
Exact frozen runtime bytes were restored; completed-output verification exited
before solver/hierarchy construction. All standard conversion and independent
reader checks then passed, with one original numerical session and hierarchy.
No inversions were repeated. The complete pipeline wall in the table includes
this downtime and is **not** an uninterrupted performance measurement; the
original numerical solve wall and its FLOP-rate denominator are unchanged.
Recovery active seconds: 320.561. Full provenance and archived failure records are
in `../cfg1020-t6-recovery-20261121/`.
