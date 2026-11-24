# Configuration 1020: standard perambulator at source/sink t=7

**Completed and verified, 2026-11-24 UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one solver session with the recorded idle windows.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | 35,399.796 s (9 h 50 min 0 s) |
| Native qprop calls, including native full FP64 checks | 33,207.556 s |
| Additional independent full FP64 checks | 1,983.192 s |
| Point extraction/write/readback | 32.237 s |
| Maximum full relative residual, either check | 2.600983065210918e-14 |
| Iterations per physical RHS | 73–82 |
| Maximum sampled solve process RSS | 5.956 GiB |
| Nominal dispatched GPU operator work | 8,338,430,984,478,720 FLOP |
| Nominal GPU operator work / total solve wall | 235.55 GFLOP/s |
| Standard projection plus complete native validation | 281.084 s |
| Complex-double GEMM time | 36.304 s |
| Maximum sampled projection process RSS | 11.278 GiB |
| Full projected gamma5 relative Frobenius discrepancy | 1.0797177139532e-14 |
| Maximum absolute gamma5 discrepancy | 2.32037276154694e-15 |
| Maximum independent QDP contraction error, norm-scaled | 7.7847080110581e-15 |
| Complete local pipeline, including whole-file hashes and reader | 82,091.898 s (22 h 48 min 12 s) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_7.sdb1020
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
7 in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is 6,907,158,528 bytes; SHA256:
`d76d9889b510f2bceba1f6118db47d109da6816823c03087f1f3957d8dcb4c6b`.

The accepted intermediate point file is `/Users/admin/Downloads/MetalClover-results/cfg1020-t7-unsmeared.npy`,
6,879,707,264 bytes, SHA256
`91f24162a20b1e6313b5d69d5bc7237e67f6ef8472a0b828eea3f2b076ff2614`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1020/eigs.t_7.mod1020`,
215,332,803 bytes, SHA256
`e4dead2322f0e76e5b5fb428771747dd268ab02e58168aa98a1d965a9464accf`. Their header reports 12³,
decay_dir=3, t_slice=7 and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to 7 in `action.xml`. The mass −0.0743,
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
and complete logs/residual evidence remain. See [the retirement record](../cfg1020-perambulator-t7-5184/intermediate-retirement.json).


## Timing and final stop

The solve-session wall and nominal FLOP-rate denominator include two recorded
idle windows totalling 5031.860 seconds. The complete pipeline wall also includes
the filesystem-access stall during postverification. These raw measurements
were preserved; they are not uninterrupted performance measurements. Both
windows resumed the same native solver process and MG hierarchy. Actual pause
and resume times are recorded in `manual-idle-window.json` and the linked
control records.

On November 24 the user requested finishing at t0=7. Final acceptance metadata
and the retained standard file's full SHA256 were rechecked successfully.
The waiting t0=8–10 workers and runtime/storage supervisor exited; their queue
entries are cancelled and no subsequent inversion started. All standard
perambulators and original inputs remain. See
`../cfg1020-stop-after-t7-20261124/` for final checksum and stop confirmation.
