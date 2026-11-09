# Checkpoint and resume record

The benchmark paused again on November 6, 2026, with **252 of 288** measured
four-RHS tiles accepted. The benchmark and its wake lock have exited. All
Metal cases and CPU multigrid and BiCGStab are complete. CPU smoother-only
FGMRES has accepted 48 RHS and would resume at offset 48. Resume only on
request. The [measured comparison](THREE_LEVEL_RESULTS.md) includes all
completed prefixes.

The benchmark and its caffeinate helper were stopped before the laptop reboot.
The run resumed on November 6, 2026, after verifying the saved source and binary
hashes. Use the progress command below for current status.

At the first pause, 195 of 288 measured four-RHS tiles were checkpointed.
All Metal cases and CPU three-level MG had finished. CPU BiCGStab had
accepted offsets 0, 4, and 8 and was interrupted at offset 12.
CPU smoother-only FGMRES was still pending. Every accepted tile passed the
full original-equation CPU FP64 residual target 1e-10.

The resume command runs from this directory and appends to the existing log:

```sh
VECLIB_MAXIMUM_THREADS=1 caffeinate -i ./build/metal_three_level \
  --size 32 --mass -.10 --csw 1 --tile-rhs 4 --rhs 96,192 --samples 1 \
  --threads 8 --coarse-iters 8 --bottom-iters 8 --resume \
  --csv results/threelevel32-compare.csv >> results/threelevel32-compare.log 2>&1
```

`python3 threelevel_progress.py` prints current progress. After all 288
measured tiles finish successfully, run
`python3 threelevel_summary.py --resumed-after-reboot` to validate all
records and traces and generate `THREE_LEVEL_RESULTS.md`.

The numerical sources/binary are frozen; hashes are in
`results/threelevel32-source-sha256.txt`. Methodology is in `THREE_LEVEL.md`.
The raw CSV and `.csv.trace.csv` must both be retained for resumption.

The measured totals retain the sleep-affected tile times. These intervals
were identified in the macOS power log on November 6 (EST).

| CPU method | Tile offset | Lid close | Full wake | Measured wall s | Krylov steps | Full FP64 residual |
|---|---:|---:|---:|---:|---:|---:|
| BiCGStab | 100 | 11:48:38 | 12:05:24 | 1346.46 | 286 | 6.91e-12 |
| Smoother FGMRES | 20 | 16:49:43 | 17:32:03 | 3565.19 | 765 | 9.53e-12 |
