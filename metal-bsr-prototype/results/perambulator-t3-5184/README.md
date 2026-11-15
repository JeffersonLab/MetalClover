# Configuration 1000: queued source/sink t=3 perambulator

All 5,184 eigenvectors and four source spins (20,736 physical RHS) are queued
for common source and sink time three. The detached local queue waits for the
entire time-two pipeline: solve, standard projection, validation and independent
DataLib reader. It also requires 17 GiB free before starting. Idle-sleep
prevention covers the queue and complete pipeline. An interrupted iPhone
connection does not stop either job.

The input is
`/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_3.mod1000`.
Its header reports 12³, decay_dir=3, t_slice=3 and 5,184 vectors. It is
215,332,763 bytes, SHA256
`b2b5b12df39c2952c3473e59ae6fb9e4decfd1b8ea6ef853b00506dbbdd1a79e`.
Stored ordering and phases are preserved; the FP32 vectors are promoted to FP64.

Only `Contractions/t_sources` changes to three in `action.xml`. FermionAction
is identical to the accepted time-zero, one and two runs: mass −0.0743,
anisotropic stout clover and antiperiodic temporal boundary conditions.
Caller/device caps remain sixteen, with one reused MG hierarchy for all RHS.
Both original full Chroma CPU FP64 residual checks must pass at 1e-12 per RHS.
No numerical code or executable changes for this run.

Output locations:

```
/Users/admin/Downloads/MetalClover-results/cfg1000-t3-unsmeared.npy
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_3.sdb1000
```

The standard output contains sixteen 5,184×5,184 complex-FP64 spin matrices.
Checks include input hashes, direct QDP contractions, complete record readback,
full gamma5 Hermiticity and independent DataLib key compatibility. The local
queue status is `cfg1000-t3-unsmeared.queue.json`; pipeline and stage status
sidecars follow the earlier runs. `staged-run.json` records inputs and settings;
`frozen-files.json` records source, executable, library and orchestration hashes.
Exact orchestration scripts are in `production/`.

With explicit user approval, the intermediate time-zero and time-one
`cfg1000-t*-unsmeared.npy` files were removed to reclaim 12.8 GiB. Their hashes
and both retained standard perambulator hashes were checked before removal.
Every standard perambulator, original gauge/eigenvector input, log, checksum
and accepted residual record is retained. See `intermediate-cleanup.json`.
Large output data remain outside Git.
