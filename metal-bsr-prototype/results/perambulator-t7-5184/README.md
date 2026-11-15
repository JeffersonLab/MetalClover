# Configuration 1000: queued source/sink t=7

All 5,184 supplied eigenvectors and four source spins (20,736 RHS) are staged
for common source and sink time 7. This job waits for complete time-6
solve, projection, verification, independent reader and queue finalization.
The action is identical to the accepted earlier runs; only the contraction
source-time entry changes.

Input: `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_7.mod1000`.
Header: 12³, decay_dir=3, t_slice=7, num_vecs=5184.
Size: 215,332,729 bytes; SHA256:
`e04dff95a03c7219613dc2feb66042a4015b61f36b501b6d54a26eaad01c76ef`.
Stored ordering/phases are preserved, promoting supplied FP32 values to FP64.

Output standard perambulator:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_7.sdb1000
```

The final file contains sixteen 5,184×5,184 complex-FP64 matrices with both
actual key times equal to 7. Caller/device caps are sixteen, one MG hierarchy
is reused across the collection, and both original full Chroma CPU FP64 checks
remain at 1e-12. Complete point checks, direct QDP contractions, all standard
record readbacks, full gamma5 Hermiticity and the independent DataLib reader
must pass before the next job starts.

`action.xml`, `staged-run.json`, `frozen-files.json` and `production/` record the
exact setup. The detached worker uses local private copies of the validated
numerical executables, libraries and controllers. It saves all completion
evidence automatically and retires only accepted point intermediates after
point/standard SHA256 rechecks. Every standard perambulator is retained.
Phone/chat disconnection does not stop the local job.

See [the five-job queue record](../perambulator-t6-t10-5184/README.md) for the
storage forecast, preservation audit, queue code and validation.
