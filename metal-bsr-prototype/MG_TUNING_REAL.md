# Metal MG tuning on the real 24³×256 light configuration

The benchmark uses the same gauge field, header action, deterministic RHS,
full-system FP64 residual definition, and `1e-12` acceptance threshold as
[the original comparison](PERFORMANCE_REAL_FIELDS.md). The final run builds
one hierarchy and reuses it for two RHS, with one RHS resident at a time.
Source 20261107 was used for tuning; source 20261108 is a holdout.
The machine is an Apple M1 Max with 32 GiB unified memory, running macOS
26.7.1. CPU residual checks use eight workers.

## Measured result

The tuned MG solver is **2.37× faster than fresh Metal BiCGStab**,
with **72.5% less fine Schur work**. Both GPU methods
solve the identical two sources from zero, including every full CPU FP64
residual check. No inversion timings overlap.

| Method | Mean solve (s/RHS) | Fine Schur applications/RHS | Maximum original Chroma FP64 full residual |
|---|---:|---:|---:|
| Previous MG (prior run) | 795.980 | 11260.5 | 2.6440e-14 |
| Fresh Metal BiCGStab | 546.759 | 8435.0 | 2.4260e-14 |
| Tuned Metal MG | 230.268 | 2318.5 | 2.3847e-14 |

The two MG solves took **231.264 s / 229.271 s**, versus **543.185 s /
550.333 s** for fresh BiCGStab. Their fine Schur counts were **2317 / 2320**,
versus **8416 / 8454**. MG used **751 / 752** FGMRES iterations; BiCGStab used
**4208 / 4227** iterations. Both methods required five defect corrections.
Relative to the earlier MG run, the tuned solver reduces fine work by
79.4% and solve time by a factor of 3.46.

The one reused hierarchy took **59.955 s**, including construction,
upload, and MG workspace setup. Adding that one cost over the two measured
RHS gives **260.245 s/RHS**, still **2.10× faster** than BiCGStab.
For N RHS, the arithmetic amortization of these measurements is
`230.267673 + 59.954768/N` seconds/RHS. Setup would contribute 0.625 s/RHS
at N=96 and 0.312 s/RHS at N=192; those are projections, not measured
96/192-RHS batch timings. Setup filtering used 1024 additional fine Schur
applications, included in setup time and excluded from solve work counts.
The discarded tuning screens are separately recorded research costs.

All four new GPU fields passed the original Chroma full-system FP64 check
at `1e-12`. The largest tuned-MG full-field relative solution difference
from the original Chroma CPU inversion was **8.6444e-12**. The standalone
FP64 full-operator probe matched the original exported Chroma probe to
1.6837e-16; the independent native EO/full probe matched to 1.7385e-16.
The native verifier also rechecked both reused CPU reference fields.

MG's full CPU FP64 refinement checks cost 10.139 s over both RHS,
2.20% of MG solve time. Their cost remains in every solve
number above. Independent original Chroma post-verification is outside
both methods' solve timings.

At MG-ready, Metal allocated 13.888 GiB; the primary CPU FP64 arrays account
for another 7.594 GiB, about 21.482 GiB combined before other host geometry,
OS, and applications. One RHS is resident. The MG buffers were released
before BiCGStab: the device allocation dropped to 3.283 GiB before its
workspaces were allocated.

[Structured summary](results/mg-tune24-final/summary.json),
[complete solve timings and work counts](results/mg-tune24-final/test.metal-solves.csv),
[full refinement traces](results/mg-tune24-final/test.trace.csv),
[independent Chroma checks](results/mg-tune24-final/test.chroma-checks.csv),
[settings and measured source commit](results/mg-tune24-final/settings.json), and
[hashed verification inputs](results/mg-tune24-final/files.sha256) preserve
the evidence. The final selected-defaults/no-screen path also passed a
12³×96 real-field smoke solve at full FP64 residual 1.4705e-14.

## Selected configuration

| Setting | Previous MG | Tuned MG |
|---|---|---|
| Levels | 3 | 3 |
| Fine aggregate | 4×4×4×8 | 3×3×3×8 |
| Fine basis vectors per chirality | 12 | 8 |
| Fine basis filtering | 8 CPU FP64 Jacobi sweeps | 128 GPU FP32 Jacobi sweeps |
| Basis filter weight | 0.65 | 0.65 |
| Second aggregate | 2×2×2×2 | 2×2×2×2 |
| Second basis vectors per chirality | 12 | 8 |
| Second basis filtering | 8 CPU FP64 sweeps | 16 CPU FP64 sweeps |
| Fine smoother | One pre + one post block Jacobi, weight 1 | Same |
| Fine FGMRES restart | 8 | 12 |
| First coarse FGMRES restart / limit / tolerance | 8 / 16 / 0.1 | 8 / 32 / 0.1 |
| Bottom BiCGStab limit / tolerance | 16 / 0.1 | 4 / 0.1 |
| First coarse dimension | 165888 | 262144 |
| Bottom dimension | 10368 | 16384 |

All block inverses and fine operators retain the original conventions.
The odd Schur operator is `A_oo - D_oe A_ee^-1 D_eo`; the clover blocks are
6×6 in each chirality. Approximate near-null vectors are filtered in FP32
on Metal; aggregate-local two-pass QR and sparse Galerkin assembly use CPU
FP64. The second basis is filtered in CPU FP64. Filtering changes only the
preconditioner; the action and full FP64 operator are unchanged.

The input is
`szscl21_24_256_b1p50_t_x4p300_um0p0856_sm0p0743_n1p265_cfg_1000a.lime`,
SHA-256 `7a4e2f5f41838bce746adf928f31b8d4f1814c904402289892fda6d20b8d9003`.
It uses mass -0.0856, xi0 4.3, nu 1.265, spatial clover coefficient
1.58932722549812, temporal coefficient 0.902783591772098, and two spatial
stout steps at rho 0.14. Time is antiperiodic. These are the header action
parameters; they were held fixed throughout tuning.

## Screening and validation

The completed screens are recorded in
[64-sweep initial tests](results/mg-tune24-a/test.screens.csv),
[64-sweep weight/bottom-limit tests](results/mg-tune24-b/test.screens.csv), and
[128-sweep tests](results/mg-tune24-final/test.screens.csv).
Each screen uses two defect corrections on RHS 0. Its partial FP64
residual and projected runtime rank candidates; they are not accepted full
solutions. One unnecessary third screen in the initial run was stopped,
with the reason preserved in its result directory.

Post-only smoothing required more Krylov iterations and more elapsed time
than one pre/post step. With 64 sweeps, weight 1.2 reduced fine work modestly.
At 128 sweeps, weight 1.0 had the lower measured runtime, while weight 1.2
used slightly fewer fine applications. The final selection follows elapsed
time. Shortening the bottom limit from 32 to 4 reproduced the same 64-sweep
iterates and residuals. Bottom solves used at most three active iterations
in the 128-sweep screens. The bottom tolerance refers to its left-preconditioned
system; its separate true FP32 residual is recorded as well.

Both complete GPU solvers use FP32 inner tolerance 0.002. Each defect
correction reconstructs the full solution, accumulates it in FP64, and
checks `||b-Mx||_2/||b||_2` over every site, color, and spin on the CPU.
Every saved final field is also checked independently by the original
Chroma FP64 operator and compared componentwise with the existing original
Chroma CPU inversion of that same source. The old CPU inversions are reused;
the Metal BiCGStab timings are fresh.

Solve time includes reduction, reconstruction, FP64 accumulation, and all
CPU FP64 residual checks. Import, common operator setup, shader compilation,
source generation, saved-field I/O, and independent native verification
are outside solve timing. Hierarchy construction/upload and its amortization
are reported separately. Discarded screening runs are tuning costs. MG
resources are released before the fresh BiCGStab runs. No inversions overlap.
Fine Schur counts include smoothing, outer products, and FP32 residual
refreshes; setup filtering is counted separately. RHS/reconstruction hopping
operations and full CPU FP64 checks are outside that Schur count but inside
solve time.

## Reproduction

`mg_real_tune` is standalone and leaves Chroma and SuperbBLAS unchanged.
The production defaults use the selected configuration and skip tuning
screens. Use `--screen-outers 2` to compare cycles; each cycle is encoded as
`pre_steps:post_steps:weight:bottom_limit`, separated by commas. The optional
fourth entry overrides `--bottom-limit` for that cycle. A zero
screen count selects the first configuration directly. The source of the
measured final run is commit `e54ead3`, identified along with executable and
source hashes in its settings; the later CLI default change only avoids
repeating discarded screens.

```sh
make -f Makefile.production build/mg_real_tune
cp build/mg_real_tune /private/tmp/mg-real-tune
mkdir -p RESULTS
cp results/performance24/action.xml RESULTS/action.xml
VECLIB_MAXIMUM_THREADS=1 python3 run_guarded.py \
  --log RESULTS/run.log --resources RESULTS/resources.csv -- \
  /private/tmp/mg-real-tune build/performance24 RESULTS/test \
  --rhs 2 --filter-sweeps 128 --cycles 1:1:1:4 \
  --screen-outers 0 --finish 1 --bicg 1
python3 verify_staged.py RESULTS --bridge build/performance24 --rhs 2 \
  --reference-directory results/performance24
```

Use a fresh output directory. Executables and verification data are staged
outside iCloud/FileProvider; hashes identify every staged file. Large binary
fields and bridges remain generated local artifacts, excluded from Git.
The solver accepts additional logical RHS sequentially without rebuilding
the hierarchy; storage and total runtime must be budgeted for the requested
number of saved fields. The completed comparison measures two RHS, rather
than claiming a measured 96/192-RHS throughput.
