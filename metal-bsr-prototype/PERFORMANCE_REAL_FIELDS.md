# Real anisotropic clover performance comparison

This document preserves the initial MG comparison. Subsequent basis,
aggregation, and smoother tuning is recorded in
[the real-field MG tuning report](MG_TUNING_REAL.md).
The subsequent [LDL/streaming implementation](LOW_MEMORY.md) fits and solves
the real 32³×256 light field; the earlier resident-layout infeasibility
calculation below is preserved as historical accounting.

The completed 24³×256 comparison used the header light mass **-0.0856** on an
**Apple M1 Max with 32 GiB unified memory**. Two deterministic dense Z4 RHS
were solved from zero, one resident at a time. All six saved fields passed the
independent original Chroma FP64 full-system check at the requested `1e-12`.

| Method | RHS 0 solve (s) | RHS 1 solve (s) | Mean solve (s/RHS) | CPU/solve speedup | Maximum native FP64 residual |
|---|---:|---:|---:|---:|---:|
| Original Chroma CPU reliable BiCGStab | 3450.684 | 3300.085 | 3375.384 | 1.00× | 4.0973e-14 |
| Metal BiCGStab | 541.207 | 542.534 | 541.870 | 6.23× | 2.4260e-14 |
| Metal three-level MG / FGMRES | 797.319 | 794.642 | 795.980 | 4.24× | 2.6440e-14 |

The single reused MG hierarchy cost **154.545 s**, including CPU construction
and device upload/workspace setup. Amortized over these two RHS, MG took
**873.253 s/RHS** and achieved **3.87×** the native CPU speed. With many RHS,
its measured solve cost approaches 795.980 s/RHS. In this configuration,
Metal BiCGStab remained faster even after amortizing the hierarchy.

MG used 3618 / 3587 FGMRES iterations versus 4208 / 4227 BiCGStab iterations,
but its fine smoother and outer iteration applied the fine Schur operator
more often: **22,521 applications** over two RHS versus **16,870** for BiCGStab.
Thus MG reduced Krylov iterations by 14.6%, increased fine Schur work by
33.5%, and increased solve time by 46.9%. This basis and smoother are not yet
competitive with GPU BiCGStab on this real gauge field. Counts cover FP32
fine Schur applications during solves; CPU FP64 filtering during hierarchy
construction is accounted for in setup time instead.

The independent native `StoutFermState` comparison had **zero componentwise
link error**, including the boundary factors. The standalone full FP64
operator probe agreed with original Chroma to **1.6837e-16** relative; the
Metal/CPU FP32 Schur probe differed by **7.5979e-8**. The largest full-field
relative solution difference from native Chroma was **1.0978e-11** for Metal
BiCGStab and **9.1697e-12** for MG. This comparison uses the original Chroma
inverter and full fields, with the same source and action for each method.

Full CPU FP64 residual checks cost 1.676 / 1.383 s for BiCGStab and
3.186 / 2.927 s for MG, less than 0.4% of the GPU solve times. The native CPU
baseline uses the installed Chroma/QDP implementation, GCC 14.2 / `-O2`,
8 OpenMP threads, and the M1 target. Swap was active during the CPU run;
these timings characterize that build on this laptop. One one-second
read-only stack sample was taken during the first CPU solve. Chroma's
reliable solver reports a zero-based `n_count`; executed iteration counts
in the CSV add one, giving 3489 / 3350 for the CPU references.

[Structured results](results/performance24/summary.json),
[native full residuals and field differences](results/performance24/test.chroma-checks.csv),
[GPU timing and work counts](results/performance24/test.metal-solves.csv),
[full refinement traces](results/performance24/test.trace.csv), and
[verified input/output hashes](results/performance24/files.sha256) preserve
the numerical evidence. Raw binary fields remain local generated artifacts.

The 32³×256 header was then read from
`szscl21_32_256_b1p50_t_x4p300_um0p0860_sm0p0743_n1p265_cfg_1000a.lime`.
It confirms **-0.0860** and the same anisotropy, clover coefficients, spatial
stout parameters, and antiperiodic time boundary. No 32³×256 solver execution,
performance number, or residual result is claimed.

| Live storage for 32³×256, one RHS, current defaults | GiB |
|---|---:|
| CPU FP64 physical gauge | 4.500 |
| CPU FP64 even/odd clover blocks | 9.000 |
| CPU FP64 source, solution, residual | 4.500 |
| Fine Metal arrays, basis and FGMRES workspace | 25.781 |
| **Lower bound before coarse operators and other costs** | **43.781** |
| Physical unified memory | 32.000 |

Even the theoretical positive minimum (restart 1 and one vector per
chirality) requires **34.406 GiB** before coarse matrices, host geometry,
setup temporaries, OS, or other apps. Its convergence was not tested; this
is a storage lower bound. The measured 24³×256 principal live storage was
19.652 GiB: 12.058 GiB in Metal buffers plus 7.594 GiB of primary CPU FP64
arrays. These measurements support the array accounting.

The device query reports a 24.960 GiB recommended GPU working set and an
18.720 GiB single-buffer maximum. The largest known default 32³×256 buffer
is the 4.500 GiB transfer basis, within the single-buffer limit; total
unified-memory capacity prevents a fully resident comparison. One RHS is
already the minimum tile. Fitting this larger case needs a smaller fixed
representation and/or streaming of FP64 checks and workspace; allocating
the current full resident layout would exceed the machine's RAM.

[32³×256 input parameters](results/performance32/input.json),
[default memory budget](results/performance32/default-memory-budget.json),
[minimum-layout bound](results/performance32/minimum-memory-budget.json), and
[feasibility decision](results/performance32/feasibility.json) record this
negative fit result without forcing an oversized allocation.

The first independent verifier launch hit macOS code-signature cache errors
after FileProvider restored its executable; the first CPU reference was
also an iCloud placeholder. Verification was completed using byte-identical
copies outside FileProvider, with all staged files hashed. The recovered
CPU reference matched the SHA-256 recorded when it was originally solved.
No solver timing was repeated or changed. `verify_staged.py` and the runner
now handle these local storage issues; staging, I/O, and verification are
outside solve timings. The failed attempts and the macOS diagnostic remain
with the results for provenance.

The performance drivers are separate from Chroma and SuperbBLAS. The native
baseline compiles Chroma's original `InvBiCGStabReliable` and FP64
`EvenOddPrecCloverLinOp`, and uses its original fixed-FP32
`EvenOddPrecDumbCloverFLinOp`. It solves the odd Schur system, performs native
even/odd reduction and reconstruction, and checks the full unpreconditioned
FP64 residual. This is an actual Chroma inverter, rather than the prototype's
CPU implementation of the same matrix.

`metal-bicgstab` is GPU BiCGStab: its odd Schur operator, Krylov vector work,
reductions, and even/odd reduction/reconstruction execute on Metal in FP32.
The full solution and defect correction are maintained on the CPU in FP64.
`metal-mg` uses GPU FGMRES with a reused three-level aggregation hierarchy.
Both GPU methods start from zero and use FP32 inner solves to relative 0.002;
full FP64 checks and corrections continue until `||b-Mx||/||b|| <= 1e-12`.
The Chroma reliable solver uses Delta 0.1 and Schur target `5e-14`, followed by
the same hard full-system `1e-12` check.

The 24³×256 input is
`szscl21_24_256_b1p50_t_x4p300_um0p0856_sm0p0743_n1p265_cfg_1000a.lime`.
The selected light mass is **-0.0856**. All action values retain the digits in
its LIME header: xi0 4.3, nu 1.265, spatial clover 1.58932722549812,
temporal clover 0.902783591772098, and two spatial stout steps at rho 0.14.
The fermion boundary is `(1,1,1,-1)`. `_per` in these ensemble names denotes
antiperiodic time. No eigenvectors or historical perambulator are used in this
timing experiment.

Two logical RHS use unit-norm dense complex Z4 sources with SplitMix64 tags
20261107 and 20261108. Each method regenerates the identical source. One RHS
is resident at a time, and the hierarchy is built once and reused. Source
generation, import, operator setup, shader compilation, saved-field I/O, and
subsequent independent Chroma verification are outside solve times. Solve
times include reduction, reconstruction, FP64 updates, and every full FP64
residual check. Hierarchy setup and its amortization are reported separately.
The large Metal run counts actual fine Schur applications, including fine
smoother and FGMRES residual work, separately from Krylov iterations. Coarse
operators and the individual hopping operations in RHS reduction and
reconstruction are outside this count.

The hierarchy uses fine aggregates 4×4×4×8, second aggregates 2×2×2×2,
12 vectors per chirality on each transfer, eight filtering sweeps with a
two-vector filter tile, one pre/post block-Jacobi step, fine FGMRES restart 8,
first-coarse FGMRES restart 8 / cap 16 / tolerance 0.1, and bottom BiCGStab
cap 16 / tolerance 0.1. These are recorded in the summary rather than inferred
from the general solver defaults.

Build against the archived source matching the installed library:

```sh
make -f Makefile.production build/performance_compare build/performance_chroma \
  build/correctness_chroma \
  CHROMA_CONFIG=/Users/admin/Downloads/macos-chroma/install/chroma/bin/chroma-config \
  CHROMA_SOURCE=/Users/admin/Downloads/qcd-macos-share-20261005/chroma.devel
```

`correctness12_prepare.py --gauge FILE --mass=-0.0856 --out-dir DIR` selects
the header action, checks agreement among matching flavors, and saves its
parameters and original XML records. `production_export FILE DIR/action.xml -
BRIDGE_PREFIX` creates the FP64 stout bridge without eigenvectors.

Run the CPU and GPU timings sequentially, with no overlapping inversions:

```sh
OMP_NUM_THREADS=8 python3 run_guarded.py \
  --log DIR/chroma.log --resources DIR/chroma-resources.csv -- \
  ./build/performance_chroma DIR/action.xml BRIDGE_PREFIX DIR/test cpu 2
VECLIB_MAXIMUM_THREADS=1 python3 run_guarded.py \
  --log DIR/metal.log --resources DIR/metal-resources.csv -- \
  ./build/performance_compare BRIDGE_PREFIX DIR/test 2 4,4,4,8
python3 verify_staged.py DIR --bridge BRIDGE_PREFIX --rhs 2
python3 performance_report.py DIR --rhs 2
```

Before accepting the large run, a native `correctness_chroma FILE
DIR/action.xml BRIDGE_PREFIX UNUSED state-only` check independently reads the
LIME field with QIO checksums, reunitarizes it, builds Chroma's original
`StoutFermState`, and compares every link including boundary conditions.
The inversion driver compares the original native full clover action with
the exported full-operator probe. The Metal driver also checks its full FP64
operator and its FP32 Schur kernel. Finally, the native verifier checks every
saved solution with the original Chroma full FP64 operator and compares all
components against the actual Chroma inversion.

The completed small validation is in `results/performance-smoke`: native CPU
18.1423 s, Metal BiCGStab 2.3861 s, and Metal MG 7.0137 s plus 6.9070 s setup.
All three full native FP64 residuals are below 4e-14; GPU solution differences
from native Chroma are below 1.2e-12 relative. These are validation numbers on
12³×96, not predictions for the larger fields.

Code, XML, logs, and numerical summaries belong in Git. Large binary gauge
bridges and saved FP64 solution fields are local generated artifacts excluded
by `.gitignore`; input hashes and saved-field hashes identify their provenance.

The recommended complete sequence is `python3 run_real_performance.py DIR
--bridge BRIDGE_PREFIX --rhs 2`. It stages byte-identical executable copies
outside iCloud, runs CPU and GPU inversions sequentially, stages and hashes
verification inputs, and requires every independent native FP64 check before
creating the summary. `verify_staged.py` can finish verification after a
storage-related failure without repeating completed solver timings.
