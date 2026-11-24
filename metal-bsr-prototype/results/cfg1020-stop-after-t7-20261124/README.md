# Completed cfg1020 time seven; stopped before time eight

The user requested finishing t0=7 production and stopping, with no t0=8 run.
The queued t8, t9 and t10 workers and runtime/storage supervisor were identified
by PID, start time and exact command, sent SIGTERM, and verified exited. Queue
entries are cancelled, persistent stop markers and automatic-start gates are
disabled, and none of these three times started a pipeline or inversion.

Time seven completed all 20,736 RHS in one numerical session/hierarchy. Both
original full Chroma CPU FP64 checks passed at 1e-12, with maximum full relative
residual 2.6009830652109176e-14. Completed-output verification, all sixteen full standard
record readbacks, direct QDP checks, gamma5 and independent DataLib reader passed.
The accepted point intermediate was retired only after proof capture and full
point/standard checksum rechecks. No additional inversion was performed.

The retained standard file is:

```
/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_7.sdb1020
```

It contains sixteen matrices of 5184×5184 complexFP64 entries, 6,907,158,528 bytes.
A fresh complete-file SHA256 check matched the accepted native output:
`d76d9889b510f2bceba1f6118db47d109da6816823c03087f1f3957d8dcb4c6b`.
The final proof is `time7-final-verification.json`. Complete numerical evidence
remains in `../cfg1020-perambulator-t7-5184/`.

The earlier filesystem stall cleared before this finalization; its underlying
cause was not established. Raw solve wall includes both recorded idle intervals;
raw pipeline wall also includes the file-access stall. Neither was rewritten.
Every standard perambulator, original input, checkpoint journal, log and durable
runtime backup remains. Source-zero cross-time blocks are still deferred.
Do not restart times 8..10 without a new explicit user request.
