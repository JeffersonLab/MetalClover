# Real anisotropic clover performance comparison

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
OMP_NUM_THREADS=8 ./build/performance_chroma \
  DIR/action.xml BRIDGE_PREFIX DIR/test verify 2 > DIR/verify.log
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
