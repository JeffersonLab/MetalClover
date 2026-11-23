# cfg1020: time-six recovery and November 23 resume

The user explicitly resumed the diagonal production queue on November 23, 2026,
after pausing it on November 21 to free disk. Disconnecting the phone or chat
does not stop these detached local processes. The scope remains P(6,6) through
P(10,10); cfg1000 times 11–20 are cancelled and cfg1020 P(t,0) is deferred.

The original time-six native solve finished all 20,736 RHS on November 20 with
exit zero, one solve session and one MG hierarchy. Both original full Chroma
CPU FP64 residual checks passed at 1e-12; the maximum full relative residual
was 2.586792573545587e-14. Its controller subsequently failed when the final
checksum verifier's shell launcher was missing. The saved numerical output,
checkpoint, per-column journal and successful native session remained intact.

On November 23, every surviving frozen runtime file matched its recorded hash.
Missing launchers, controllers, source snapshots and per-time staging files
were restored byte-for-byte from their recorded preimages/workspace archives.
No executable, library, numerical source, action parameter or solver setting
changed. `restoration-20261123.json` records all restored and backed-up files.
An independent durable copy is stored outside temporary storage at:

```
/Users/admin/Downloads/MetalClover-runtime/cfg1020-frozen-20261123
```

`resume_20261123.py recover6` invokes the original native command with `--resume`
on the **already complete** checkpoint. This checks every saved column and
exits before constructing a solver or hierarchy. The original one-session
summary is then produced, followed by the original standard projection,
complete record readbacks, direct QDP contractions, gamma5 check and independent
DataLib reader. Original failed statuses and logs are preserved under
`time6-before-recovery-20261123/`. The successful recovery record is
`time6-recovered.json`; it is written only after full proof capture and accepted
point retirement. No inversions are repeated.

The numerical solve wall remains the original 30,090.162 seconds. The recovered
pipeline's elapsed wall includes the intervening runtime failure and user pause;
it is explicitly marked and must not be interpreted as an uninterrupted
performance measurement. Recovery active time is recorded separately.

Time-seven through time-ten workers call the unchanged original batch driver.
Their only scheduling change is a separate storage forecast through their own
time. This permits a fully funded prefix to run while later jobs wait for disk;
the complete-sequence forecast still truthfully reports any deficit. Original
18 GiB start, 4 GiB reserve and native projection guards remain. All standard
perambulators and original inputs are retained. Only accepted point intermediates
are retired after proof capture, full SHA256 rechecks and an open-file check.

At resume, about 30.3 GiB was free after the durable backup. The forecasts fit
through time eight; times nine and ten required approximately 2.4 and 9.2 GiB
additional free space, respectively. These are live capacity estimates, not
completion promises. A runtime supervisor refreshes temporary-file timestamps
and restores only missing files from verified durable copies. A present file
whose hash differs is rejected instead of being overwritten.

`resume-launch-20261123.json` records detached process identities and log paths.
`runtime-supervisor.json` is the live runtime heartbeat. Original pause records
are retained as history. Completed scientific evidence is stored in the normal
`../cfg1020-perambulator-tN-5184/` records. The three recovery fixtures test
exact restoration, rejection of corrupt/present-mismatched files, and permitting
only the safe storage prefix; the original nine queue-safety and four allocation
tests also passed. These control fixtures do not replace the production checks.

Large output payloads and the durable runtime backup remain outside Git.


Recovery completed successfully on November 23: both original full FP64 checks,
every saved point-column checksum, all sixteen standard record readbacks, direct
QDP contractions and the independent DataLib reader passed. Full gamma5 relative
Frobenius discrepancy was 1.08440794058118e-14. The accepted point intermediate
was retired only after the final full SHA256 rechecks; the standard remains.
See `time6-recovered.json` and the normal time-six validation record.
