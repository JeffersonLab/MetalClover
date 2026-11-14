# Configuration 1000: queued source/sink t=2 perambulator

The same full 5,184-eigenvector, four-source-spin calculation is staged for
source and sink time two (20,736 physical RHS). It starts only after the
time-one solve, standard projection, complete verification and independent
analysis-reader key check finish successfully. The queue and both pipelines
run locally on the Mac, with idle-sleep prevention, independently of the
user's iPhone connection. An interrupted phone connection is not a stop request.

The selected file is
`/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_2.mod1000`.
Its header reports 12³, t_slice=2, decay_dir=3 and 5,184 vectors. The file is
215,332,741 bytes; SHA256 is
`2155ec707251ab97f6c1145aa43f3f7178ae66bd10d52bdac5681541470d6b71`.
Ordering and phases are preserved, with FP32 sources promoted to FP64.
`action.xml` differs from the time-zero action only in `Contractions/t_sources`.
The mass −0.0743, anisotropies, clover coefficients, spatial stout action and
antiperiodic temporal boundary conditions are identical.

Production uses caller/device caps of sixteen, one hierarchy reused across
all 20,736 RHS, and both original Chroma CPU FP64 residual checks at 1e-12.
The final standard perambulator is sixteen 5,184×5,184 complex-FP64 matrices.
The projection checks input hashes, direct QDP contractions, complete record
readback, full gamma5 Hermiticity and independent DataLib reader compatibility.
No solver code or numerical settings change between times one and two.

Output locations:

```
/Users/admin/Downloads/MetalClover-results/cfg1000-t2-unsmeared.npy
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_2.sdb1000
```

Live queue state is `cfg1000-t2-unsmeared.queue.json`; pipeline state is
`cfg1000-t2-unsmeared.pipeline.json` once execution starts. Native progress is
in the point-prefix `.run.json`. Each stage produces a `.complete.json` only
after its checks pass. `staged-run.json` records inputs and settings;
`frozen-files.json` records tested source/executable/library hashes. The queue
rechecks the frozen files and original inputs before starting time two.

`production/queue_after.py` and `production/full_run.py` are the exact local
job orchestration scripts. They wait for the accepted time-one result and
then run the existing solve, projection and independent-reader controllers.
`temporary-artifact-cleanup.json` records obsolete derived benchmark data
reclaimed only after checking archived hashes and confirming the files were
not open. All original gauge/eigenvector inputs and requested production
outputs are retained. Large output data remain outside Git.
