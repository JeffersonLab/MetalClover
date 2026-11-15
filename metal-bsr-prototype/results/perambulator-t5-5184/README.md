# Configuration 1000: queued source/sink t=5 perambulator

All 5,184 supplied eigenvectors and four source spins (20,736 physical RHS)
are staged for common source and sink time five. The detached local queue
waits for successful completion of both the time-four pipeline and its queue:
solve, standard projection, complete verification and independent DataLib reader.
Idle-sleep prevention covers the wait and all stages. An interrupted connection,
including an iPhone disconnect, does not stop the local jobs.

Input eigenvectors:
`/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_5.mod1000`.
The header reports 12³, t_slice=5, decay_dir=3 and 5,184 vectors. The file is
215,332,870 bytes, SHA256
`8c73450b08f42806cb3db973a35714a7dd566fbe914d82182ef34a324a18a2e0`.
Stored ordering and phases are preserved, promoting the supplied FP32 values
to FP64.

Only `Contractions/t_sources` changes to five in `action.xml`. All fermion-action
parameters remain identical: mass −0.0743, anisotropic clover, two spatial stout
steps with rho=0.14 and orthog_dir=3, and antiperiodic temporal boundary conditions.
The same frozen native binaries and library are used, with caller/device caps
of sixteen, one MG hierarchy reused across all RHS, and both original full
Chroma CPU FP64 residual checks at 1e-12. No numerical solver code changes.

Output locations:

```
/Users/admin/Downloads/MetalClover-results/cfg1000-t5-unsmeared.npy
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_5.sdb1000
```

The final standard perambulator contains sixteen 5,184×5,184 complex-FP64
matrices. Projection verifies input hashes, direct QDP contractions, complete
record readback, full gamma5 Hermiticity and independent DataLib key compatibility.
Queue status is `cfg1000-t5-unsmeared.queue.json`; pipeline and stage sidecars
follow the earlier runs. `staged-run.json` records inputs and settings;
`frozen-files.json` records source, binary, library and orchestration hashes.
Exact local scripts are in `production/`.

Every standard perambulator and original input is retained. Four obsolete
generated temporary benchmark payloads were reclaimed only after matching
their archived SHA256 hashes and confirming they were not open; this does not
affect the running time-four job. See `temporary-artifact-cleanup.json`.
The queue saves all accepted time-four evidence into its workspace run record
before retiring that run's unsmeared `.npy` intermediate if space is needed.
Retirement requires complete acceptance, point/standard SHA256 rechecks and
confirmation that the intermediate is not open; it creates a separate
`intermediate-cleanup.json`. The queue requires at least 17 GiB free before
starting time five. No unfinished point file is removed. Large data stay
outside Git.
