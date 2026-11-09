# Clover-Wilson odd Schur complement and Metal FP32 elimination

This extension adds a gauge-derived clover term with **c_sw=1.0**, uses Chroma's
odd-site Schur convention, and moves RHS reduction and even-site reconstruction
onto Metal in FP32. The full solution, accumulation, and true-residual checks
remain CPU FP64. See [CLOVER_RESULTS.md](CLOVER_RESULTS.md) for measurements at
96 and 192 RHS on the M1 Max. Chroma, QDP, and SuperbBLAS checkouts remain
read-only dependencies for inspection, not linked libraries.

The clover operator and CPU/Metal Schur products were implemented and validated
**before** changing the elimination precision. The first-stage validation log is
`results/clover-stage1.log`; its driver is retained as `clover_stage1.mm`.

## Exact Chroma/QDP convention

The following local sources determine the convention:

* `/Users/admin/Downloads/chroma.devel/lib/meas/glue/mesfield.cc`: four oriented
  plaquettes and an anti-Hermitian field with factor 0.125.
* `/Users/admin/Downloads/chroma.devel/lib/actions/ferm/linop/clover_term_qdp_w.h`:
  coefficient scaling by 0.5 and `makeClovSiteLoop`'s two packed chiral blocks.
* `/Users/admin/Downloads/chroma.devel/lib/actions/ferm/linop/eoprec_clover_linop_w.cc`:
  even-block inversion and the retained odd Schur operator.
* `/Users/admin/Downloads/qdpxx.devel/include/qdp_primspinvec.h`:
  the actual Gamma(1), Gamma(2), Gamma(4), Gamma(8), and Gamma(15) actions.

The new driver uses those QDP gamma matrices, with
`gamma5=diag(1,1,-1,-1)`. They are the complex conjugates of the explicitly
documented basis used by the earlier Wilson prototype. Only the clover driver
selects this QDP basis; historical drivers retain their existing conventions.
Spin basis and clover coefficient are therefore consistent with each other.

For each site and mu<nu, the four positively oriented leaves are

```text
Q_mu_nu(x) = loop(+mu,+nu,-mu,-nu)
           + loop(+nu,-mu,-nu,+mu)
           + loop(-mu,-nu,+mu,+nu)
           + loop(-nu,+mu,+nu,-mu).

F_AH_mu_nu = (Q_mu_nu - Q_mu_nu^H) / 8
C(x)       = -(c_sw/2) sum_(mu<nu) gamma_mu gamma_nu F_AH_mu_nu(x)
A(x)       = (mass+4) I + C(x).
```

`F_AH` is anti-Hermitian, not the Hermitian physical F used in some texts.
No trace subtraction is added, matching the inspected `mesField` implementation.
Gauge loops are periodic; antiperiodic fermion boundary factors apply to hopping,
not to this gauge-field construction.

The **1/4 versus 1/2** depends on whether both ordered index pairs or only mu<nu
are summed and on the definition of sigma. With
`sigma_mu_nu = [gamma_mu,gamma_nu]/2`, the same explicit convention is
`C = -(c_sw/4) sum_(all mu,nu) sigma_mu_nu F_AH_mu_nu`.
For only mu<nu it is `-(c_sw/2) sum sigma F_AH`. Definitions that move factors
of i between sigma and F change the displayed sign/factors. The implementation
avoids that ambiguity by checking its gamma-product construction against an
independent transcription of Chroma's packed-component formulas.

This is the isotropic mass-form operator with c_sw=1 and u0=1. There is no hidden
kappa rescaling. The input positive bare masses are 0.1 and 0.5 in the benchmark.

## Two 6×6 matrices on every site

C commutes with diagonal gamma5, so each local A has two Hermitian 6×6 blocks:

```text
A(x) = A_plus(x) (+) A_minus(x)
```

This applies at both even and odd sites. Within a block, the six indices are
`q = 3*spin_within_chirality + color`, matching Chroma's packed spin/color order.
The external field layout stays color × four spins × RHS with RHS fastest.
Full 6×6 matrices are stored in this prototype; triangular packing is not used.

The FP64 A_ee blocks are inverted using complex Cholesky during setup, separately
for each chirality and even site. A nonpositive pivot causes an error. The
inverses are reused for all RHS and iterations. FP32 copies of A_oo and A_ee^-1
are rounded from these authoritative FP64 matrices. Runtime FP32 application
does not imply FP32 or GPU factorization.

## Chroma's odd-site Schur form

Write the full operator in checkerboard order as

```text
M = [ A_ee  H_eo ]
    [ H_oe  A_oo ],       H = -0.5 Dslash.
```

The retained system and reconstruction are exactly

```text
S_o       = A_oo - H_oe A_ee^-1 H_eo
g_o       = b_o - H_oe A_ee^-1 b_e
S_o x_o   = g_o
x_e       = A_ee^-1 (b_e - H_eo x_o).
```

There is **no scalar mass division or extra Schur normalization**. If D denotes
raw Dslash rather than the off-diagonal matrix block H, the product has the
equivalent factor 1/4 from the two -1/2 hopping coefficients, as in Chroma.

CPU and GPU Schur products perform an odd-to-even hop, local A_ee^-1, an
even-to-odd hop, and `A_oo*x_o - product`. The local operations are genuine
chiral 6×6 products. Packed even and odd maps require even periodic extents.
The clover term does not change the off-diagonal hopping parity.

The outer mixed-precision loop applies this construction to the **current full
FP64 residual** to obtain a correction. It accumulates both reconstructed
parities into an FP64 solution and accepts the result only when every RHS meets
`||b-M64*x64||/||b|| <= 1e-10`. Zero RHS use an absolute residual criterion.
This is a full-system check, not a recursive or reduced-system residual.

## Five measured paths

* `cpu_fp64`: FP64 Schur BiCGStab, FP64 RHS reduction/reconstruction.
* `cpu_mixed_host64`: FP32 CPU Schur solve, FP64 CPU elimination.
* `cpu_mixed_fp32`: FP32 CPU Schur solve and FP32 CPU elimination.
* `metal_mixed_host64`: FP32 Metal Schur solve, FP64 CPU elimination. This
  preserves the earlier transfer policy while adding the correct clover operator.
* `metal_mixed_fp32`: FP32 Metal Schur solve, reduction, and reconstruction.

All use the identical clover operator and RHS within each case. Mixed inner
tolerance is 0.002. FP64 inner tolerance is 1e-10. All full residuals and
accumulated solutions are FP64 on the CPU. Complex BiCGStab acts directly on
the nonsymmetric Schur operator, with no normal equations. Inner iteration
limit is 512 and outer correction limit is 12. Failures and nonfinite results
produce an error, not a successful timing row.

The new device path performs:

1. CPU FP64-to-FP32 conversion of the full residual into a shared Metal buffer.
2. GPU gather of even/odd sites, A_ee^-1 application, H_oe, and RHS subtraction.
3. GPU-resident BiCGStab initialization directly from the reduced RHS. The reduced
   RHS is never downloaded or promoted to FP64 to initialize the inner solver.
4. GPU H_eo, subtraction from b_e/current r_e, A_ee^-1 application, and scatter of
   the full FP32 correction, followed by completion.
5. CPU FP64 accumulation and full FP64 clover-Wilson residual evaluation.

The CPU cast in step 1 remains necessary because Metal has no native FP64 data
arithmetic. Shared memory avoids a separate output download, but synchronization
at CPU/GPU boundaries is explicit and timed. Four BiCGStab iterations are encoded
per command buffer by default, including padded work after convergence.

## Build, validate, and benchmark

Run from this directory with Xcode tools and GPU access:

```sh
make
make clover-test
make clover-bench
python3 clover_summary.py results/clover.csv --output CLOVER_RESULTS.md
```

The benchmark defaults are equivalent to

```sh
./build/metal_clover --sizes 8 --rhs 96,192 --mass .1,.5 --csw 1 \
  --gauge-strength .25 --samples 3 --threads 8 --chunk 4 \
  --tol 1e-10 --inner-tol .002 --csv results/clover.csv
```

A quick complete comparison:

```sh
./build/metal_clover --sizes 4 --rhs 3,7 --mass .25 --samples 1 \
  --csv results/clover-quick.csv
```

Use `--periodic-time` to change the fermion time boundary, `--csw 0` to remove
clover, and `--gauge-strength 0` for a free field. The original first-stage
driver can be rebuilt with `make build/clover_stage1` and run separately.

## Timing and validation

Warm total time includes resets, input conversion, RHS preparation, all Schur
BiCGStab operations/reductions, command encoding/submission/waits, reconstruction,
FP64 updates, and full CPU FP64 residual checks. Setup includes field/block
construction and inversion, allocation/upload, and manufactured data; it is
recorded separately. Each timed solve restarts from zero, reuses prepared
operators/workspaces, and uses the same RHS. There is one full warmup and three
trials per method, with rotating order. CPU paths request eight GCD workers.

`inner_ms` includes `prepare_ms` and `reconstruct_ms`; subtract these to isolate
the inner solve, which also includes its initialization. The report does this
without double counting. Counts are summed maximum per-RHS iterations across
outer steps; GPU dispatched counts also include command-chunk padding.

Each benchmark RHS has a seeded known FP64 solution and b generated with the
original FP64 operator. Known-solution error is checked after timing, with
relative error limited to 1e-7. The old Wilson mass-only condition bound is not
assumed for clover. Every component of every true residual is checked inside
timing, with no sampling.

Validation includes an analytic quantized Abelian field: four equal leaves give
`F_01 = diag(i*sin(theta), -i*sin(theta), 0)`. This independently fixes the field
orientation and 1/8 normalization. Other checks cover Chroma component formulas,
chiral separation, Hermiticity, independent Gaussian-versus-Cholesky inverses,
local gauge covariance, direct full operator and odd Schur products, CPU/Metal
agreement, and direct device RHS/reconstruction against known solutions.

The final eight-case self-test exercises all five solver paths at small even
volumes with odd RHS, a zero RHS, free field, and c_sw=0. The staged validation
log demonstrates clover correctness before the device elimination change.
Original Wilson/even-odd regression drivers remain available.

## Reuse and scope

Clover field strengths are constructed once from the gauge links. Local clover
matrices and their even-site inverses are then reused across iterations and
many RHS. This supplies additional local matrix reuse and arithmetic intensity;
it does not automatically reduce the hopping kernel's gauge-link traffic or
guarantee that adding clover makes a solve faster.

The matrices still come from synthetic near-identity SU(3) links, not a production
ensemble. This test is isotropic, positive bare mass, and unpreconditioned beyond
even/odd; no multigrid or tuned half-spinor Dslash is included. Neither the CPU
nor GPU implementation is claimed to be optimal. FP64 residual work can still
limit total GPU speedup, particularly as RHS and working sets grow.

Files: `clover_common.hpp` (field, blocks, inverses, operators and elimination),
`clover.metal` (local block and transfer kernels), `clover_checks.hpp` (independent
checks), `clover.mm` (comparison driver), and `clover_summary.py` (report).
