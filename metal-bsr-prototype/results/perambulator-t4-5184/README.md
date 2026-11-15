# Configuration 1000: running source/sink t=4 perambulator

**Started automatically on 2026-11-15 at 17:29:49 UTC**, after complete
time-three validation, evidence capture and verified intermediate cleanup.
The initial source batches have passed both full CPU FP64 residual checks;
`handoff-validation.json` records the accepted handoff, not a completed run.

All 5,184 supplied eigenvectors and four source spins (20,736 physical RHS)
are running for common source and sink time four, using the same procedure as
the earlier runs. The detached local queue waits for the complete time-three
pipeline: solve, standard projection, complete verification and independent
DataLib reader. Sleep prevention covers the queue and all stages. Connection
loss, including an interrupted iPhone connection, does not stop the local jobs.

Input eigenvectors:
`/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_4.mod1000`.
The header reports 12³, t_slice=4, decay_dir=3 and 5,184 vectors. The file is
215,332,828 bytes, SHA256
`5605ffb6568cc8fd3800c66cf09bbe21be740abdf254612f73624261372a6114`.
Stored ordering and phases are preserved, promoting the supplied FP32 values
to FP64.

Only `Contractions/t_sources` changes to four in `action.xml`. The fermion
action is identical: mass −0.0743, anisotropic stout clover and antiperiodic
temporal boundary conditions. The same frozen native binaries and library
are used, with caller/device caps of sixteen, one MG hierarchy reused across
all RHS, and both original full Chroma CPU FP64 residual checks at 1e-12.
No numerical code, executable or solver setting changes.

Output locations:

```
/Users/admin/Downloads/MetalClover-results/cfg1000-t4-unsmeared.npy
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_4.sdb1000
```

The final standard perambulator contains sixteen 5,184×5,184 complex-FP64
matrices. Projection verifies input hashes, direct QDP contractions, complete
record readback, full gamma5 Hermiticity and independent DataLib key
compatibility. Queue status is `cfg1000-t4-unsmeared.queue.json`; pipeline and
stage sidecars follow the earlier runs. `staged-run.json` records all inputs;
`frozen-files.json` records source, binary, library and orchestration hashes.

Following the user's request to repeat the previous procedure, all standard
perambulators and original inputs are retained. If space is needed, only the
accepted time-two and time-three unsmeared `.npy` intermediates are reclaimed.
This occurs after time-three completion, after its complete evidence is saved
in the repository workspace, after rechecking both intermediate and retained
standard-file SHA256 hashes, and after confirming neither intermediate is
open. The queue then requires at least 17 GiB free before starting time four.
The actual cleanup produces `intermediate-cleanup.json`; it has not happened
at staging. No point file for an unfinished solve is removed.

Exact local scripts are in `production/`. The dependency evidence helpers
capture all accepted time-three logs and checks before its intermediate can
be retired; the large data files remain outside Git.
