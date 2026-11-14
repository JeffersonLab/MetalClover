# Standard perambulator from the accepted point-sink solutions

`project_perambulator.cc` forms

```
tau[spin_sink, spin_source](sink_vector, source_vector)
    = sum_(x,colour) conj(V[x,colour,sink_vector])
                        * solution[source_vector,spin_source,x,colour,spin_sink]
```

at source and sink time zero. Both ends use the supplied MOD eigenvector
ordering and phases. Eigenvectors stored as complex FP32 are promoted to FP64;
the contraction and stored values are complex FP64. The operation uses CPU
Accelerate complex-double GEMM and requires no further inversions or hierarchy
construction. The verified unsmeared collection remains its input.

The output is Chroma's usual FILEDB/SDB `propElemOp` representation, using the
same `KeyPropElementalOperator_t`, `ValPropElementalOperator_t` and record-store
writer as the restored distillation measurement. For this single source/sink
time and mass label `light`, it has sixteen spin-block records. Each value is
a 5,184 by 5,184 matrix with sink-vector rows and source-vector columns.
Keys contain `t_source=0`, `t_slice=0`, `spin_src`, `spin_snk`, and `mass_label`.
The metadata carries the anisotropic stout-clover action and the original
point-file/solver provenance. Spatial phasing is zero.

The original full four-dimensional CPU FP64 residual checks at 1e-12 are
inherited from the accepted solutions. This postprocessing verifies input
SHA256 hashes before reading those solutions, performs 128 independent QDP
sink contractions, decodes all sixteen complete output records with Chroma's
value reader and checks bitwise serialization roundtrips, then checks the
full projected matrix's gamma5 Hermiticity. Only after those checks does it
rename the temporary database to its final filename. The controller records
the output's complete-file SHA256.

## Build and run

With a native Chroma installation:

```sh
make -C metal-bsr-prototype build/project_perambulator \
  CHROMA_CONFIG=/path/to/chroma-config
python3 metal-bsr-prototype/run_projected.py \
  --binary /path/to/project_perambulator \
  --point-prefix /path/to/cfg1000-t0-unsmeared \
  --eigenvectors /path/to/eigs.t_0.mod1000 \
  --output /path/to/prop.n5184.light.t0_0.sdb1000 --vectors 5184
```

As with the other native helpers, a build-tree configuration may need
`CHROMA_EXTRA_CPPFLAGS` for its Chroma and dependent-library include paths.
The existing validated build can instead compile this source through its
native test-program flags and libraries. The converter is separate from the
MetalClover library and Eloy's measurement; it does not require changes to
either numerical implementation.

The controller uses Python's standard library, with `OMP_NUM_THREADS=1` and
`VECLIB_MAXIMUM_THREADS=8` for CPU BLAS. It prevents idle sleep, samples RSS
and free disk, and stops above 26 GiB process RSS or below 4 GiB free disk.
Choose a fresh output filename. An interrupted conversion leaves `.partial`
and diagnostic sidecars; it is not a resumable inversion run. Its accepted
source solutions are unchanged and can be projected again.

`OUTPUT.projection.csv` records each spin block's GEMM time, QDP contraction
error and encoded-value checksum. `OUTPUT.projection.xml` contains the native
validation measurements and database metadata. `OUTPUT.complete.json` records
successful completion, source/output hashes and inherited full residual checks.

The independent analysis-side compatibility check uses
`Datalib::TensorStorage::StorageProp4`. On the four-vector pilot it reads all
sixteen spin blocks into the standard Prop4 tensor, and every complex element
matches the Chroma writer's value bytes exactly. See
[the run record](results/standard-peram-t0-5184/README.md) for the full result.


## Completed full-basis file

The configuration-1000 standard perambulator is complete at:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_0.sdb1000
```

It contains sixteen 5,184×5,184 complex-FP64 matrices (6.433 GiB), with both
times equal to zero. Projection plus full record and symmetry validation took
271.028 s, including 35.262 s of complex-double GEMM. The largest direct QDP
contraction error was 5.7913e-15 (normalized by the two vector norms), and the
full matrix's gamma5 Hermiticity discrepancy was 1.0800e-14 in relative
Frobenius norm. All full matrix records passed bitwise readback checks.
Input SHA256 checks confirm that the original accepted 20,736 solutions and
their supplied eigenvectors were used. Their full Chroma CPU FP64 residual
checks remain documented in the original solve record.
