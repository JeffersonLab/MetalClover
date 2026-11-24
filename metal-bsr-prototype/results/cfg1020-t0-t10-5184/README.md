# Configuration 1020: diagonal perambulators at times 1–10

**Queued and launch-verified on 2026-11-17.** The user selected diagonal blocks
first: P(1,1), P(2,2), …, P(10,10). The cfg1000 time-eleven through time-twenty
queue and its storage supervisor were cancelled before any of those inversions
started. The active cfg1000 calculation and its time-ten successor are unchanged.
Cfg1020 time one waits for complete cfg1000 time-ten solve, projection, evidence
capture, full validation and independent reader completion; subsequent cfg1020
times follow sequentially. These are queued jobs, not completed numerical results.

The input gauge is
`/Users/admin/Downloads/unitary_distillation/cfgs/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1020.lime`.
All eight fermion-action entries in its LIME header agree with the validated
cfg1000 action: mass −0.0743, clover coefficients 1.58954576550983 and
0.902866347885228, xi0=4.3, nu=1.265, two spatial stout steps with rho=0.14,
orthog_dir=3, and fermion boundaries 1 1 1 −1. The lattice is 12³×96, with
FP32 gauge input promoted through the unchanged native Chroma path. Complete
SciDAC payload checksums match suma=0a2a46cb, sumb=1ad9a52c. `gauge-metadata.json`
records the full input hash, LIME record inventory and action comparison.

Every supplied eigenvector file for source/sink times zero through ten has
**5,184 finite eigenvalues in `MODMetaData/Weights`**, stored in descending
order. For t=0 they span 0.168540958176172 to 11.8314590418238. Stored vector
ordering, phases and weights are preserved; no eigenvalue sorting is applied.
Input hashes, per-slice extrema and metadata are in `eigenvector-inputs.json`.
The existing native exporter embeds the complete eigenvector XML in its point
metadata; standard conversion retains that metadata, including the weights.

Each requested diagonal block uses all 5,184 eigenvectors and four source
spins: 20,736 physical RHS, with caller/device caps sixteen and one MG hierarchy
reused for the whole collection. Both original full four-dimensional Chroma
CPU FP64 residual checks remain at 1e-12 per RHS. Standard projection retains
direct QDP contraction checks, all sixteen complete record readbacks, full
gamma5 Hermiticity and the independent DataLib key reader. Complete evidence
is saved automatically, including the last job, before an accepted point
intermediate is retired after point/standard SHA256 rechecks and an open-file
check. Every standard perambulator and original input is retained.

The numerical executables, Metal library, dependent libraries and controllers
are unchanged. They use the frozen private runtime at
`/private/tmp/chroma-times6-10-production/runtime`; keep it through the entire
cfg1020 sequence. The new per-time orchestration is under
`/private/tmp/chroma-cfg1020-timeN-production`, with the driver/storage monitor in
`/private/tmp/chroma-cfg1020-production`. No rebuild or new GPU solve was needed
for staging, and no cfg1020 GPU work starts before the cfg1000 dependency.
Nine synthetic queue/retirement tests, four allocation tests, actual native
argument-list checks, Python syntax checks and final frozen-file hash checks
passed. These are control checks, not cfg1020 convergence results.

The storage forecast passed at launch with about 92.5 GiB free. It includes
finishing cfg1000 times nine/ten, keeping every standard output, retaining all
ten cfg1020 diagonal blocks, retiring only accepted point intermediates,
1 GiB of global reserve, 256 MiB of evidence/logs per remaining job, standard
metadata/padding, an 18 GiB inversion-start gate and the native projection's
space requirement. It predicted about 10.8 GiB free at the tightest point and
17.3 GiB after final point retirement. The live monitor checks every thirty
seconds; no file deletion/movement is performed by that monitor. See
`disk-budget-diagonal-at-launch.json` for the complete snapshot and assumptions.

Final standard files remain under `/Users/admin/Downloads/MetalClover-results`,
ending in `.prop.n5184.light.t0_N.sdb1020`. Complete per-time workspace records
are `results/cfg1020-perambulator-tN-5184/`. Queue/driver PIDs are recorded in
`queue-launch.json`. `python3 batch_status.py` reads current states. Detached
workers and idle-sleep prevention span waiting and all numerical stages;
phone/chat disconnection does not stop them. The laptop must remain running.
Future completed evidence is saved locally but is not automatically committed.

## Source-zero blocks: checked, not yet computed

No existing cfg1020 perambulator was found in the checked output/input paths.
Earlier cfg1000 jobs retained **only P(t,t)**. In particular, their source-zero
output retained sink time zero; it does not contain P(1,0) through P(10,0).
Those off-diagonal blocks cannot be recovered from the saved diagonal blocks.

Computing all ten P(t,0) blocks needs one source-time-zero bank of 20,736 RHS,
with additional sink extraction/projections rather than extra inversions for
each sink. The current frozen exporter supports equal source/sink times, so
a multi-sink export extension is needed before that bank runs. It can preserve
the full 4D FP64 residual checks. A cross-time gamma5 check must compare P(t,0)
with P(0,t); self-Hermiticity of P(t,0) alone is not the appropriate check.

Ten additional standard blocks occupy approximately 64.3 GiB; keeping those
and the ten diagonal blocks requires approximately 128.7 GiB in final files
alone. A conservative future pipeline with ten point spools and one standard
conversion at a time needs about 58 GiB more free space than the current
post-diagonal forecast. This is a capacity/algorithm feasibility estimate,
not a numerical test of an implemented cross-time exporter. The user chose
to defer these source-zero blocks. See `source0-feasibility.json`.

## November 23 recovery and resume

Cfg1000 through time ten and cfg1020 diagonal times one through five completed
with their full acceptance checks. Time six finished all 20,736 numerical
inversions, but a missing temporary launcher stopped its subsequent verifier.
The user paused the queue on November 21, then explicitly resumed it on
November 23. Exact frozen runtime files were restored and backed up outside
temporary storage. Time six resumes verification and standard conversion
without any new inversion session or hierarchy; times seven through ten follow
as detached local jobs. All original numerical settings and checks remain.

Each remaining worker now reads a scoped storage forecast through its own time,
allowing a fully funded prefix to proceed while later jobs wait. The separate
full-sequence forecast continues to report the actual remaining deficit. At
resume the forecasts fit through time eight, with about 9.2 GiB additional free
space needed for the complete sequence through ten. No standard output or
original input is removed. See
[the recovery record](../cfg1020-t6-recovery-20261121/README.md) for exact hashes,
original failure evidence, detached process records and the durable runtime
backup. Phone/chat disconnection does not stop the local jobs.


## Final state: stopped after time seven, November 24

Cfg1020 diagonal times 1..7 completed all production checks. On the user's
explicit instruction, the waiting times 8..10 were cancelled before any pipeline
or inversion started; their workers and runtime/storage supervisor exited.
All standard perambulators and inputs are retained. The final time-seven standard
file passed a fresh complete SHA256 check. No further automatic production is
authorized. See [final verification and stop proof](../cfg1020-stop-after-t7-20261124/README.md).
