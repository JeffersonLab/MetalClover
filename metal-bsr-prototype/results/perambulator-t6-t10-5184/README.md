# Configuration 1000: sequential source/sink times 6–10

**Queued and launch-verified on 2026-11-15 UTC.** Five full
standard-perambulator jobs are queued after time five, in order
6 → 7 → 8 → 9 → 10. Each uses all 5,184 supplied eigenvectors and four source
spins: 20,736 physical RHS per time, 103,680 additional RHS in total. Source
and sink times are equal in each job. These are staged jobs, not completed
numerical results.

The action and numerical executables are unchanged from the accepted earlier
times: mass −0.0743, anisotropic clover, two spatial stout steps with rho=0.14,
orthog_dir=3 and antiperiodic temporal boundary conditions. Only
`Contractions/t_sources` changes in each action XML. Caller/device RHS caps
remain sixteen, with one MG hierarchy reused for the complete collection at
each time. Both original full four-dimensional Chroma CPU FP64 residual
checks remain at 1e-12 for every RHS.

Every job waits for both the preceding pipeline and its queue finalization.
Acceptance includes complete point-column checksum verification, projection
into sixteen 5,184×5,184 complex-FP64 matrices, direct QDP contraction checks,
all sixteen full record readbacks, full gamma5 Hermiticity and the independent
DataLib key reader. A failure stops advancement; there is no automatic restart
that could silently build another hierarchy. The final time-ten job also saves
its complete evidence without requiring another job or an interactive session.

The native programs, Metal library, shared-library dependencies and two existing
controllers have private, hash-verified copies under
`/private/tmp/chroma-times6-10-production/runtime`. The native binaries and
libraries are byte-identical to the validated originals. Small shell wrappers
select the private libraries. CPU-only usage probes verified actual dynamic
loading; no additional GPU solve was launched for those probes. Per-time
controllers, actions and manifests live under
`/private/tmp/chroma-timeN-production`. Later workspace edits do not change the
frozen numerical runtime for this queue.

The workers are detached local processes with idle-sleep prevention. Phone or
chat disconnection does not stop them. The laptop must remain running; the
queue does not recover automatically across shutdowns or reboots. Read current
states with `python3 batch_status.py`, or inspect the queue/pipeline sidecars in
`/Users/admin/Downloads/MetalClover-results`.

All standard perambulators, original inputs and complete validation evidence
are retained. Only accepted unsmeared `.npy` intermediates are retired, after
evidence capture, whole-file point/standard SHA256 rechecks and an open-file
check. Each retirement is recorded separately. An 18 GiB start gate and a
forecast covering every remaining output protect the standard projection's
space requirement. `disk-budget.json` records the forecast and its reserves.

Additional storage was recovered through transparent APFS compression of old,
closed generated diagnostics and build artifacts. Their paths and complete
byte contents are preserved, with SHA256 verification before and after atomic
replacement. Original inputs, numerical datasets, source files, current runtime
and perambulators were excluded. See `storage-preservation.json`, the per-file
`storage-preservation.jsonl` and tool provenance. Storage preparation overlapped
the time-four calculation, so that wall time includes concurrent CPU/I/O work.
The supplemental `binary-preservation/` audit covers old executables with
incompressible individual chunks; at least 10% net savings and unchanged full
SHA256 are still required. `binary-compression-probe.json` records that probe.

`queue-safety-validation.json` covers synthetic dependency/storage and retirement
guard fixtures; it makes no numerical solver claim. `runtime-validation.json`+covers argument lists, frozen hashes and library loading. Individual run records
contain eigenvector headers/hashes, action XML and scripts:

- [Time six](../perambulator-t6-5184/README.md)
- [Time seven](../perambulator-t7-5184/README.md)
- [Time eight](../perambulator-t8-5184/README.md)
- [Time nine](../perambulator-t9-5184/README.md)
- [Time ten](../perambulator-t10-5184/README.md)

Final files stay outside Git, under `/Users/admin/Downloads/MetalClover-results`,
with names ending in `.prop.n5184.light.t0_N.sdb1000`. Future completed evidence
is written to the corresponding workspace run directory; the detached workers
do not perform automatic Git commits.

The final preparation audit passed with 55.035 GiB actually free. The forecast
reserves 1 GiB globally and 256 MiB for logs/evidence at each remaining time,
plus up to 64 MiB of metadata/padding per standard file. It predicts at least
20.124 GiB before any new inversion and 6.997 GiB at the tightest output peak,
with all standards retained. About 13.404 GiB remains after the final accepted
point intermediate is retired. See `storage-audit-validation.json`,
`disk-budget.json`, `queue-launch.json` and `queue-status-after-launch.json`.
The final forecast assumes no additional large unrelated disk writes.
