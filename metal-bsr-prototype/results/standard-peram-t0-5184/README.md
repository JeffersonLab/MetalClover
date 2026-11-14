# Configuration 1000: standard time-zero perambulator

This is the standard `V_dagger M_inverse V` perambulator for all 5,184 supplied
color eigenvectors, with source and sink both at t=0 on the real 12³×96 field.
The mass label is `light`; the quark mass is −0.0743, with the same anisotropic
stout-clover action and antiperiodic temporal boundary condition as the accepted
unsmeared solve. Stored eigenvector ordering and phases are preserved.

**Completed and independently verified on 2026-11-14.** No additional
inversions or hierarchy builds were needed.

| Measurement | Result |
| --- | ---: |
| Source / sink eigenvectors | 5,184 / 5,184 |
| Spin records | 16 |
| Matrix values | complex FP64 |
| Native projection plus complete validation | 271.028 s |
| CPU complex-double GEMM time | 35.262 s |
| Maximum sampled process RSS (5 s sampling) | 11.310 GiB |
| Largest independent QDP contraction error, norm-scaled | 5.7913e-15 |
| Full matrix gamma5 Hermiticity relative Frobenius discrepancy | 1.0800e-14 |
| Maximum absolute gamma5 discrepancy | 1.5019e-15 |
| Database bytes | 6,907,158,528 (6.433 GiB) |

The native time includes basis loading, projection, standard serialization,
complete Chroma record readback, bitwise value checks and full gamma5 symmetry
checking. Source-file SHA256 preflight and final whole-database SHA256 are
outside that native time. The sixteen GEMMs perform a nominal 8*n*n*k*16
real FLOPs, with n=k=5,184; those matrix multiplications account for only part
of the conversion/verification time. BLAS used eight threads, with OpenMP set
to one.

The output database SHA256 is
`d987f2a56d393e8941d5315c8d842762f22c6299b2cd506947cb060c5acabef4`.
`full.complete.json`, `full.projection.xml`, `full.projection.csv`, `full.log`
and `full.resources.csv` record the complete conversion and checks.
`production-validation.json` records both input integrity and final acceptance.
The independent DataLib consumer also recognizes the full-size file's one
Prop4 key at source/sink t=0 with mass label `light`; all matrix values were
read back by Chroma, while the separate small consumer test checks complete
Prop4 spin/vector mapping.

The standard Chroma FILEDB/SDB file is retained at:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_0.sdb1000
```

There are sixteen `KeyPropElementalOperator_t` / `ValPropElementalOperator_t`
records. Each matrix has 5,184 sink-vector rows and 5,184 source-vector columns;
the key identifies source/sink time, source/sink spin and mass label.
Both the matrix values and the contraction arithmetic are complex FP64.

The projection reuses the accepted 20,736 unsmeared solutions from
[the full-basis solve](../unsmeared-t0-5184/README.md). Its original full Chroma
CPU FP64 residual checks at 1e-12 are retained, with maximum relative residual
2.4063810813519627e-14. No inversion is repeated and no new hierarchy is built.
The source NPY SHA256 is checked against that completion record before conversion.

The native converter checks 128 independently evaluated QDP sink inner products,
reads and decodes every complete matrix through Chroma, checks all serialized
value bytes, and evaluates gamma5 Hermiticity over the entire projected matrix.
The output is published after those checks and its full SHA256 is recorded.

`validation/reader-validation.json` records the independent DataLib
`StorageProp4` compatibility test: all sixteen spin blocks and all 256 complex
values of a four-vector pilot match the Chroma writer's encoded values exactly.
The pilot also passes the same direct-contraction, full readback and gamma5 tests.
`validation/frozen-files.json` identifies the tested converter and executable.

See [the converter guide](../../STANDARD_PERAMBULATOR.md) for build/run commands,
formats and timing scope. Large output databases remain outside Git.
