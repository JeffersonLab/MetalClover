# Temporary idle window for New Wick contraction

Completed on 2026-11-17 Eastern. The user requested up to thirty minutes of
idle CPU for the existing benchmark. Cfg1020 time one had not started; the
active workload was cfg1000 time ten. The native solve, its resource controller
and parents, the ten queued cfg1020 start gates and the storage monitor were
kernel-suspended with SIGSTOP after complete process-identity checks. The same
processes resumed with SIGCONT; no solver exited or hierarchy was rebuilt.

The pause lasted **1544.714 seconds**, from
22:50:10 to 23:15:54 Eastern. A detached
watchdog provided a thirty-minute fallback and released the window early after
the peer's completion marker. All fifteen process resumptions were verified.
The native PID remained 46030, solver session zero, one near-null/SVD hierarchy.
The paused checkpoint was 17,104 RHS; post-resume validation observed 17,152 RHS
and a maximum full FP64 relative residual of 2.36246e-14. Cfg1020 time one
remained queued behind complete cfg1000 time-ten validation; the storage guard
passed. All output files and original numerical checks were preserved.

New Wick contraction completed all six paired benchmark runs with byte-identical
output checks. A rejected attempt caused by background CPU activity remains in
its own records. The peer retained numeric idle thresholds and documented its
recognition of verified kernel-stopped processes as inactive. Spotlight/display
handling belonged to that chat; this controller signalled only its fifteen
explicitly listed queue/solver processes. The final peer results are preserved
in `peer-results.snapshot.json`, with its SHA256 in `summary.json`.

**Timing interpretation:** native session/qprop/total wall timers may include
this pause. The original raw timing reports are not rewritten. Use the recorded
pause duration when distinguishing production wall time from active work. The
same note is saved as `manual-idle-window.json` in the cfg1000 time-ten record.
Future perambulator generation and full FP64 validation continue normally.

`plan.json`, `pause-status.json`, `controller-launch.json`, `controller.log`,
`paused-checkpoint.json`, `peer-completion.json` and `resume-validation.json`
preserve the process identities, deadline and actual resume proof. The saved
controller is historical, one-use execution evidence; do not rerun it against
later processes or import it as a library.
