# Full eigenvector-source to unsmeared-point-sink output

The separate Chroma test program `t_metal_clover_unsmeared` uses Eloy's
existing `SB::doInversion` callback without changing his projected perambulator
measurement. It solves all four source spins for each supplied color
eigenvector at a selected source time, retains unsmeared site/color/spin sinks
at that same time, and streams accepted columns to disk. The current harness
targets one rank on the real 12³×96 configuration 1000. `--time T` selects the
common source/sink time (default zero); the action XML's `t_sources` and the
supplied MOD eigenvector file must agree with it. Different source/sink times
or other geometries require an explicit extension of this harness.

For 5,184 eigenvectors this means 20,736 physical RHS, not 5,184 spinor RHS.
The solver and its MG hierarchy are constructed once and reused throughout
the run. The source-spin loop surrounds the eigenvector loop, as in Eloy's
caller. Caller `max_rhs=16` gives 1,296 array calls; `max_rhs=1` remains valid.
The device `RHSBatch` is fixed to sixteen independently of the caller cap.

Every full four-dimensional solution is checked against the original Chroma
unpreconditioned clover operator in CPU FP64, both by the native qprop and by
an independent evaluation in the output callback. Both relative residuals
must be at most 1e-12 before the selected time slice is retained. The other 95
sink times are discarded only after these checks.

## Inputs and conventions

- Gauge: `szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime`.
- Action XML: `results/unsmeared-t0-5184/action.xml`; mass −0.0743,
  spatial/temporal clover coefficients 1.58954576550983 / 0.902866347885228,
  xi0=4.3, nu=1.265, stout rho=0.14 with two spatial steps,
  and antiperiodic temporal boundary conditions.
- Sources: `/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_0.mod1000`.
  Its 5,184 stored complex-FP32 eigenvectors are promoted to FP64 without
  sorting or phase changes. Their stored eigenvalue order is descending.
  The inversion tolerance is relative to these supplied sources; it does
  not increase the precision of the input eigenvectors.

The gauge/action state is constructed by Chroma, including stout smearing and
anisotropic clover conventions. The source basis may use different smearing
from the fermion action. This output has point sinks and eigenvector sources;
it is not the sink-eigenvector-projected standard perambulator.

## Build and run

Configure the Chroma `rge/devel-datalib-reorg` branch with SuperbBLAS and
`--with-metal=/path/to/metal-prefix`, then explicitly build:

```sh
make -C mainprogs/tests t_metal_clover_unsmeared
python3 run_unsmeared.py \
  --binary /path/to/chroma-build/mainprogs/tests/t_metal_clover_unsmeared \
  --gauge /path/to/configuration.lime --action /path/to/action.xml \
  --eigenvectors /path/to/eigs.t_0.mod1000 \
  --prefix /path/to/results/cfg1000-t0-unsmeared \
  --first 0 --vectors 5184 --max-rhs 16 --method MetalMGClover
```

The controller is Python standard-library only. It records binary/input SHA256
hashes, runs with `OMP_NUM_THREADS=1` and `VECLIB_MAXIMUM_THREADS=1`, prevents
idle sleep while solving, and samples process RSS/free disk every five seconds.
Its default guard requests a checkpointed stop above 26 GiB process RSS or
below 4 GiB free disk. The native program reserves the complete output payload
on disk before starting the solves.

Use a fresh prefix. To pause, send SIGTERM to the controller or native process;
the native writer finishes its current batch and exits 75. Repeat the exact
command with `--resume` to continue. The caller cap can change when resuming;
immutable action, source, solver and layout metadata must match. A restart
constructs a fresh hierarchy, so one uninterrupted run is required for a
single hierarchy over the entire collection.

## Output and reading

`PREFIX.npy` is little-endian complex128, C order, with shape
`(5184, 4, 12, 12, 12, 3, 4)` and axes
`(source_vector, source_spin, z, y, x, sink_colour, sink_spin)`.
The payload plus NPY header occupies 6,879,707,264 bytes (about 6.41 GiB).
The native QDP site order is converted explicitly to natural spatial order.

```python
import numpy as np
p = np.load('/path/to/results/cfg1000-t0-unsmeared.npy', mmap_mode='r')
value = p[vector, source_spin, z, y, x, colour, sink_spin]
# Matrix with point-sink rows and eigenvector/spin source columns:
matrix = p.reshape(20736, 20736).T
```

`PREFIX.meta.xml` contains source/action metadata, ordering and input hashes.
`PREFIX.columns.csv` records both full residuals, source norm, iteration count
and payload FNV1a64 checksum for each accepted column. `PREFIX.state` is the
atomic batch checkpoint. Payload/journal data are synced before publishing it.
Resume truncates uncommitted journal tails and verifies every committed
column checksum before continuing.

`PREFIX.run.json` is live status, with per-session logs and resource CSVs.
Only status `complete` together with `PREFIX.complete.json` means final
verification passed: the native completed-output resume path rereads every
column checksum, and Python checks the NPY header/size, all journal indices
and residuals, then hashes the entire file with SHA256. The completion report
includes total runtime and nominal dispatched GPU operator GFLOP/s; this
numerator excludes CPU FP64 and Krylov vector work, while the denominator
includes setup, source input, residual checks and output. Final checksum
rereading and SHA256 are outside that solve-session time.

Large eigenvectors and generated propagators stay outside Git. Only code,
action XML, compact validation evidence, hashes and reproduction instructions
are versioned.

## Completed full-basis run

The complete configuration-1000 run finished on 2026-11-14: all 20,736 RHS,
one reused hierarchy, and both full original Chroma CPU FP64 residual checks
passed at 1e-12. The largest relative residual was 2.4064e-14. Total solve-session
wall time was 8 h 22 min 20 s, with peak process RSS 5.955 GiB on the 32 GiB
M1 Max. Nominal dispatched operator work divided by this total wall time was
277.49 GFLOP/s, with the counting scope described above.

The verified payload is
`/Users/admin/Downloads/MetalClover-results/cfg1000-t0-unsmeared.npy`.
See [the full run record](results/unsmeared-t0-5184/README.md) for exact input,
output and binary hashes, all column residuals, timing details and native
completed-output checksum verification.


## Selecting another source time

For example, use `run_unsmeared.py --time 1` with
`results/perambulator-t1-5184/action.xml` and the supplied
`eigs_mod.n5184/1000/eigs.t_1.mod1000` file to solve source and sink time one.
The native program accepts the same `--time 1` option. Time zero remains the
default and existing time-zero checkpoints retain the same metadata layout.
The solver checks that the action source-time list and eigenvector metadata
agree with the selected time before creating output. The projector infers
the time from the accepted point-file metadata and writes that time into
both standard perambulator keys and XML.
