# Standard perambulator from the accepted point-sink solutions

`project_perambulator.cc` forms

```
tau[spin_sink, spin_source](sink_vector, source_vector)
    = sum_(x,colour) conj(V[x,colour,sink_vector])
                        * solution[source_vector,spin_source,x,colour,spin_sink]
```

at the common source/sink time recorded in the accepted point file. Both ends
use the supplied MOD eigenvector ordering and phases. Eigenvectors stored as
complex FP32 are promoted to FP64;
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


The full source/sink time-one repeat is also complete. All 20,736 RHS passed
both original Chroma CPU FP64 checks at 1e-12 (maximum 2.38185e-14), using
one reused hierarchy. The verified standard file ends in
`.prop.n5184.light.t0_1.sdb1000`, with source/sink time one in its actual keys.
The solve took 30,095.368 s and standard projection/validation 273.205 s;
the full projected gamma5 relative Frobenius discrepancy was 1.08110e-14.
See [the time-one run record](results/perambulator-t1-5184/README.md).
The complete time-two repeat also passed all 20,736 RHS with one hierarchy:
maximum full CPU FP64 residual 2.44056e-14, solve wall
30,198.007 s, and full projected gamma5 relative discrepancy
1.08482e-14. Its verified standard filename ends in
`.prop.n5184.light.t0_2.sdb1000`; see
[the time-two record](results/perambulator-t2-5184/README.md).

The complete time-three repeat passed all 20,736 RHS with one hierarchy:
maximum full CPU FP64 residual 2.46445e-14, solve wall
30,172.532 s, and full projected gamma5 relative discrepancy
1.08171e-14. Its verified standard filename ends in
`.prop.n5184.light.t0_3.sdb1000`; see
[the time-three record](results/perambulator-t3-5184/README.md). With explicit
user approval, the time-zero and time-one unsmeared `.npy` intermediate payloads
were retired after final standard-file validation and checksum checks to make
space. Their standard perambulators and all original inputs and logs remain.

Time four started automatically after complete time-three verification, with the same
20,736-RHS procedure and full CPU FP64 checks. The local queue also saves all
time-three evidence before reclaiming validated point intermediates when
needed for space, retaining every standard perambulator and original input.
See [the time-four record](results/perambulator-t4-5184/README.md).

The accepted time-two/time-three unsmeared `.npy` intermediates were subsequently
retired after complete evidence capture and point/standard checksum rechecks to
make room for time four. Every standard perambulator and original input remains;
see [the cleanup record](results/perambulator-t4-5184/intermediate-cleanup.json).

Time five is queued after complete time-four pipeline and reader verification,
using all 20,736 RHS, the unchanged action, one reused MG hierarchy and both
full Chroma CPU FP64 residual checks at 1e-12. The detached local queue saves
accepted time-four evidence before verified intermediate retirement when needed
for space; all standard perambulators and original inputs remain. See
[the time-five record](results/perambulator-t5-5184/README.md).


Source/sink times six through ten are queued sequentially after time five,
with 20,736 RHS at each time, the unchanged action and both full Chroma CPU
FP64 checks at 1e-12. Each job uses one reused hierarchy, saves complete evidence
and retains its standard perambulator. The detached workers continue without
a phone/chat connection. See [the five-job queue record](results/perambulator-t6-t10-5184/README.md)
for the frozen private runtime, full-sequence storage budget and checks.


Times eleven through twenty are queued sequentially after time ten, with the
same 20,736-RHS action and checks at every time. Their separate detached workers
retain every standard output and continue waiting without a phone connection.
The full-sequence storage guard currently needs approximately 65 GiB additional
free space before the extension can start; existing runs remain unchanged.
See [the extension record](results/perambulator-t11-t20-5184/README.md).


The cfg1000 time-eleven through time-twenty extension was cancelled on
2026-11-17 before those inversions began. Cfg1000 continues through time ten,
then the replacement cfg1020 queue computes P(1,1) through P(10,10) with the
same action/checks and all standard files retained. Cfg1020 eigenvalues are
present in `MODMetaData/Weights`; source-zero P(t,0) blocks are not present in
older diagonal outputs and remain deferred by the user's priority selection.
See [the replacement queue and feasibility record](results/cfg1020-t0-t10-5184/README.md).
