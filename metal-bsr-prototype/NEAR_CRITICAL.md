# Near-critical clover-Wilson comparison

This follow-up keeps the same 16⁴ synthetic SU(3) field, gauge-strength parameter
0.25, antiperiodic time, c_sw=1, QDP gamma conventions, and full CPU FP64
residual acceptance. The selected test mass is **m₀=-0.12** in the prototype's
`(m₀+4)I - Dslash/2 + C` normalization.

## Selecting the mass

Wilson fermions have a bare mass parameter; a value of zero is not a calibrated
critical mass for this synthetic background. We used a finite-configuration
spectral diagnostic to locate a low-gap region before timing solvers.
For the full operator define `Q(m)=gamma5*M(m)`. Gamma5 Hermiticity makes Q
Hermitian and `Q²=M†M`, so the smallest absolute Q eigenvalue equals the
smallest singular value of M. This gap convention follows
[Del Debbio et al., section 2](https://arxiv.org/abs/hep-lat/0512021).

`critical_probe` applies the full matrix in CPU FP64, using `M(m)=M(0)+mI`,
and runs Lanczos on Q² with full two-pass reorthogonalization. Small real
tridiagonal eigenproblems use Accelerate LAPACK. Every reported candidate is
reconstructed on the full lattice and its actual FP64 Q² eigen-residual is
computed. Both random and smooth starts were used, with a different-seed
500-step check near the observed minimum. The smooth start spans both signs
of the lowest temporal Fourier wave, plus local random noise.

The estimates put the low-gap region near bare mass -0.14. The selected -0.12
lies about 0.02 above that finite-volume proxy, with an estimated gap close to
the minimum. These are approximate Rayleigh estimates. The finite-volume
Lanczos runs did not attain the requested 1e-5 eigen-residual, and their actual
residuals and convergence flags are reported rather than claiming a precise
critical mass. A physical PCAC critical mass requires different measurements.

Antiperiodic time remains unchanged. The gap of this small near-unit-link
configuration stays finite in the low-gap region. It is therefore a
solver stress test near one fixed field's estimated gap minimum, with finite-volume
and synthetic-gauge limits. Free-field checks independently reproduce

```
sigma_min(m)^2 = min_p [(m + sum_mu(1-cos p_mu))^2 + sum_mu sin² p_mu],
p_space = 2*pi*n/L,       p_time = (2*n+1)*pi/L.
```

The existing known solution is used for RHS generation and post-solve checks;
hierarchy setup and the spectral search are independent of it.

## Making the coarse correction suitable

At m=-0.12 the earlier four-iteration coarse BiCGStab cap produced corrections
whose actual coarse residual could grow by a large factor. Fine FGMRES still
converged, but spent substantially more fine Schur work. A small-RHS sweep
increased the cap to 8, 12, and 16 while retaining the preconditioned coarse
target of 0.1. Twelve iterations achieved a good actual coarse residual on
that pilot; sixteen gave the same convergence and provides room for harder
RHS in the 96-wide batch.

The final configuration keeps twelve local basis vectors per chirality,
4⁴ aggregates, eight setup-filtering sweeps at weight 0.65, one pre/post
clover-block smoothing step at weight 1.0, and fine FGMRES restart 8. The
coarse cap is **16**, and this configuration change is part of the result.
All fine inner targets remain 0.002; the full FP64 target remains 1e-10.
BiCGStab remains the same odd-Schur baseline, with a GPU convergence check
every iteration. Smoother-only FGMRES retains the same two-step smoother and
restart. No normal equations are introduced in the three timed solvers;
Q² is used only by the separate spectral diagnostic.

The sparse matrix is still the full Galerkin `P† S_o P`, with 6,144 unknowns
and 8,448 dense 24×24 sparse blocks. The actual negative-mass 16⁴ hierarchy
passes Galerkin/gamma5 identities and strict CPU/Metal coarse-solve checks.
The local even clover factors remain positive in this configuration. Both
processors pass the negative-mass small-lattice full-solve self-test, including
zero RHS and even reconstruction.

## Timing and RHS protocol

96 RHS are resident simultaneously. Each trial solves global RHS 0–95 and
then 96–191, starting each tile from zero. The first measured tile supplies
the 96-RHS result; the actual sum of both measured tiles supplies the 192-RHS
result. These are paired prefix measurements. CPU and Metal receive identical
operator and input data. Large per-method workspaces are allocated sequentially;
the hierarchy and FP64 outer buffers are shared for the mass/volume case.

One resident-tile warmup precedes each method's measured trials. This run
starts with Metal, followed by CPU. Within each processor methods run in
blocks, retaining their workspace across trials. Ordinary local system
activity is present; no thermal or scheduling isolation is imposed.

Timing includes resets, FP32 even/odd reduction/reconstruction, all operator
and vector/reduction work, CPU small QR, GPU submissions/waits, FP64 solution
updates, and full CPU FP64 residual checks on every component/RHS. Setup and
manufactured RHS generation have separate columns outside solve time.
Known-solution error is checked after residual acceptance outside timing.
Every RHS must satisfy `||b-M64*x64||/||b|| <= 1e-10`. The spectral locator's
approximate eigen-residual does not change this solve criterion.

Screening used 12 resident RHS and numerical/operator-work criteria. Screening
times are excluded from the final timing comparison; some convergence checks
and spectral diagnostics overlapped. Final timed methods run sequentially.

## Reproduction

```sh
make build/critical_probe build/metal_solver_compare
VECLIB_MAXIMUM_THREADS=1 ./build/critical_probe --self-test --smooth-start
VECLIB_MAXIMUM_THREADS=1 ./build/critical_probe \
  --size 16 --steps 300 --tol 1e-5 --smooth-start \
  --mass .1,-.10,-.12,-.14 --csv results/critical-refine.csv
VECLIB_MAXIMUM_THREADS=1 ./build/critical_probe \
  --size 16 --steps 500 --tol 1e-5 --seed 92871 --smooth-start \
  --mass -.14 --csv results/critical-independent.csv
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \
  --self-test --mass -.12 --coarse-iters 16 \
  --csv results/nearcritical-final-self-test.csv
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \
  --sizes 16 --rhs 3 --tile-rhs 3 --mass -.12 --coarse-iters 16 --validate-only
VECLIB_MAXIMUM_THREADS=1 ./build/metal_solver_compare \
  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass -.12 --samples 3 \
  --threads 8 --csw 1 --bi-chunk 1 --coarse-iters 16 --coarse-tol .1 \
  --metal-first --csv results/nearcritical-compare16.csv
python3 nearcritical_summary.py
```

See [NEAR_CRITICAL_RESULTS.md](NEAR_CRITICAL_RESULTS.md) for measured results
and the retained spectral residuals. The earlier positive-mass comparison is
[SOLVER_COMPARISON_RESULTS.md](SOLVER_COMPARISON_RESULTS.md). Its historical
times are not used as fresh near-critical timing samples. All code remains
inside the standalone prototype; Chroma and SuperbBLAS sources are unchanged.
