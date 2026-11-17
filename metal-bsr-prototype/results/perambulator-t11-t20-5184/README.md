**Cancelled on 2026-11-17: superseded by cfg1020 P(1,1) through P(10,10).**
All ten waiting workers and their storage supervisor exited; cfg1000 continues
through time ten. No standards or inputs were removed. See
[the replacement queue](../cfg1020-t0-t10-5184/README.md).

# Configuration 1000: queued source/sink times 11–20

**Queued and launch-verified on 2026-11-16 Eastern time. Storage is pending.**
Ten detached workers wait in order 11 → 12 → … → 20 after complete time-ten
validation/finalization. Each uses all 5,184 eigenvectors and four spins:
20,736 RHS per time, 207,360 additional physical RHS. Both source and sink equal
the selected time. These are queued jobs, not completed numerical results.

The action is identical to the accepted earlier runs: mass −0.0743, anisotropic
clover, two spatial stout steps with rho=0.14, orthog_dir=3 and antiperiodic
temporal boundary conditions. Only the contraction source-time entry changes.
Caller/device caps remain sixteen. One hierarchy is reused for all RHS at
each time, with both original full four-dimensional Chroma CPU FP64 residual
checks at 1e-12 for every RHS.

The queue retains complete point-column checksums, direct QDP contractions,
all sixteen standard record readbacks, full gamma5 Hermiticity and the independent
DataLib key reader. Complete evidence is saved automatically, including time
twenty, before accepted point intermediates can be retired. Retirement rechecks
point/standard SHA256 and confirms that the point file is closed. Every standard
perambulator and original input is retained.

The new workers have separate scripts/manifests under
`/private/tmp/chroma-timeN-production`, with the queue driver and storage
supervisor in `/private/tmp/chroma-times11-20-production`. They reuse the
unchanged, hash-verified private numerical runtime under
`/private/tmp/chroma-times6-10-production/runtime`. Keep that shared runtime
until time twenty has completed. No existing worker, numerical source, library,
executable or solver setting was changed, and no new GPU test ran during staging.

All ten eigenvector headers and full hashes passed. The native argument lists,
frozen manifests and Python syntax passed; eight synthetic queue/retirement
tests and four synthetic full-sequence storage forecasts passed. These control
tests do not establish future numerical convergence. Per-time records contain
the exact actions, inputs, hashes and scripts.

## Storage hold

At staging, 34.62 GiB was free. The conservative forecast for every remaining
standard output through time twenty requires approximately 99.95 GiB free in
that initial state: **65.33 GiB additional space**. It includes a 1 GiB global
reserve, 256 MiB of logs/evidence per remaining time, standard metadata/padding,
the native projection's space check, an 18 GiB inversion start threshold and a
4 GiB minimum output reserve. See `disk-budget-initial.json`.

No external volume was mounted. A read-only inventory of old research/build
artifacts did not establish enough safe local savings for this deficit. No
additional files were deleted, moved or compressed during this extension.
The existing time-seven solve and queued times eight through ten continue
unchanged. The extension cannot start while the full-sequence storage guard
reports insufficient space. Additional storage or separately approved cleanup
is needed to release it.

The detached storage supervisor recomputes the forecast every thirty seconds,
accounting for current allocations and retaining all standard outputs. A stable
allocation snapshot is required. If enough local space becomes available, the
guard passes automatically; each worker still waits for its predecessor's
complete validation and finalization. Actual free space is checked again before
every inversion. Failed or paused predecessors do not release subsequent jobs.
The supervisor never deletes or moves files.

Current storage state:
`/private/tmp/chroma-times11-20-production/storage-supervisor.json`.
Current job states: `python3 batch_status.py` or the corresponding local copy.
`queue-launch.json` records all ten verified worker/sleep-preventer PIDs.

All workers and the supervisor run locally, detached from the phone/chat
connection. Idle-sleep prevention covers the wait and subsequent stages. The
laptop must stay running; shutdown/reboot requires manual recovery. Final
standard names end in `.prop.n5184.light.t0_N.sdb1000` under
`/Users/admin/Downloads/MetalClover-results`. Future completion evidence is saved
in each workspace run directory; the workers do not commit results automatically.
