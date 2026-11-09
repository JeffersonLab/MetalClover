# Even/odd-preconditioned Wilson BiCGStab

This adds a packed half-volume Schur-complement solve to the standalone Wilson
prototype. [EVENODD_RESULTS.md](EVENODD_RESULTS.md) compares the preconditioned
and unpreconditioned CPU FP64, CPU mixed, and Metal mixed paths at 96–192 RHS.
Chroma and SuperbBLAS are not linked or modified.

## Operator and convention

The full operator, synthetic SU(3) links, spin basis, boundary conditions, and
manufactured RHS are unchanged from [WILSON.md](WILSON.md). With `d = mass+4`
and H denoting the off-site `-0.5 Dslash` blocks, the checkerboard form is

```text
         [ d I    H_eo ]
D_full = [             ] .
         [ H_oe   d I  ]
```

Eliminating odd sites and normalizing by d gives

```text
S     = I - H_eo H_oe / d^2
g_e   = (b_e - H_eo b_o / d) / d
S x_e = g_e
x_o   = (b_o - H_oe x_e) / d .
```

The scalar-diagonal case follows the structure in `getEvenOddPrec` near line
2278 of `/Users/admin/Downloads/chroma.devel/lib/util/ferm/mgproton.cc`. Its
`use_Aee_prec` variant also has `I-H_eo H_oe/d²`; here the reduced unknown is
explicitly x_e and the RHS is scaled accordingly. A clover diagonal would require
nontrivial local block inverses and is not implemented by this scalar formula.

Even sites have even coordinate sum, not merely even linear index. Both parity
sets are packed into consecutive half-volume arrays. Each extracted hopping
row contains the 16 identity/gamma terms that connect only to the opposite
checkerboard. The diagonal term is excluded. All neighbor indices are remapped
to compact opposite-parity indices, including antiperiodic time boundary signs.

All periodic lattice extents must be even: an odd wraparound extent has
same-parity nearest neighbors. The CLI rejects odd extents before running.
The independent unpreconditioned driver still supports odd extents.

## CPU and Metal implementation

A CPU Schur application performs H_oe into a half-volume scratch vector, H_eo
back into the even output, and `y = x - y/d²`. ARM NEON kernels and eight GCD
worker tasks are used in the measured CPU paths. FP64 and FP32 versions use their
respective coefficient arrays and arithmetic.

The Metal path performs the same two half-volume hopping kernels followed by a
combination kernel, with barriers between dependent dispatches. The vectors and
scratch remain GPU resident during the inner BiCGStab solve. Complex reductions,
recurrence coefficients, and per-RHS convergence flags remain on Metal as before.
Four BiCGStab iterations per command buffer is the default. Padded iterations
and the second operator application after early convergence remain part of timing.

Both operators use the existing generic color/spin kernels, not a specialized
half-spinor Dslash. Two half-volume hopping passes are approximately one full
volume of hopping work, not half: iteration-count reduction alone is not a
prediction of total speedup.

The shared BiCGStab implementations accept an operator callback for the Schur
application; their direct full-Wilson application remains the default. This
keeps the recurrence implementation identical between the two operator choices.
`wilson_common.hpp` now shares the Wilson constructor and independent stencil
checks between the drivers. `build/metal_wilson` remains independently runnable.

## Mixed precision and full-system stopping

The outer refinement loop retains a full-volume FP64 solution and residual.
For every correction, starting with `r = b-D64*x`:

1. Gather the even and odd residuals and form the reduced RHS in **CPU FP64**:
   `g_e = (r_e - H_eo64*r_o/d)/d`.
2. Solve the half-volume Schur system with complex BiCGStab. CPU mixed and Metal
   mixed use FP32 corrections and inner tolerance 0.002; CPU FP64 uses 1e-10.
3. Promote the even correction and reconstruct the odd correction in **CPU FP64**:
   `delta_o = (r_o - H_oe64*delta_e)/d`.
4. Accumulate the full FP64 correction, apply the original full Wilson operator
   in CPU FP64, and check **every component of every RHS**.

The solve is accepted only when every RHS satisfies the original equation's
`||b-D64*x|| / ||b|| <= 1e-10`, not when a Schur or recursive FP32 residual alone
is small. Zero RHS use the same absolute residual convention as the earlier
solver. Full-system-converged RHS receive zero correction RHS in subsequent
outer steps. All-FP64 solves also undergo this full reconstruction/check process
and refine again if necessary.

Using FP64 RHS reduction and reconstruction is an explicit implementation choice
to retain accurate elimination outside the FP32 inner solve. Both mixed paths
pay those CPU costs. They can limit the net Metal gain after inner iterations
are reduced. Full FP64 residual checks remain unchanged and are inside timing.

The FP64 Wilson coefficients are the numerical authority. FP32 hops and the
rounded inverse-d-squared factor approximate that operator only inside the inner
solve. Known-solution comparisons are independent validation after each solve.
Iteration limits, breakdown, nonfinite values, missed true tolerances, or failed
known-solution checks cause a failing exit code.

## Reproduction and timing

Run from this directory with Xcode build tools and Metal access:

```sh
make
make evenodd-test
make evenodd-bench
python3 evenodd_summary.py results/evenodd.csv --output EVENODD_RESULTS.md
```

The complete comparison defaults to:

```sh
./build/metal_evenodd --sizes 8 --rhs 96,192 --mass .1,.5 \
  --samples 3 --threads 8 --chunk 4 --preconditioner both \
  --csv results/evenodd.csv
```

A quick check:

```sh
./build/metal_evenodd --sizes 4 --rhs 3,7 --mass .25 --samples 1 \
  --csv results/evenodd-quick.csv
```

Other options: `--preconditioner none|even_odd|both`, `--gauge-strength .25`,
`--tol 1e-10`, `--inner-tol .002`, and `--periodic-time`. The standalone file
loader expects to find `kernels.metal`, `bicgstab.metal`, and `evenodd.metal`
in the current directory. The inherited memory guard limits a full vector to
18 million complex elements. This comparison holds workspaces for both methods;
it is not a measurement of the minimum memory needed by either method alone.

For every case, all six methods use exactly the same full Wilson operator and
RHS. There is a full warmup per method followed by three measured solves from
zero, with rotating method order. Results are local measurements with ordinary
system activity; per-trial timings and spread are retained in the CSV.

Total solve time includes reset/copies, casts, Schur products, all BiCGStab
recurrences/reductions, Metal encoding/submission/waits, FP64 RHS reduction,
FP64 reconstruction, FP64 solution updates, and full FP64 residual checks.
Fixed RHS norms, operator/RHS generation, initial checkerboard extraction,
allocation/upload, and compilation are outside warm solve time. Setup is
reported separately. The final known-solution comparison is outside solve time.

In the raw CSV, `inner_ms` includes `eo_prepare_ms` and `eo_reconstruct_ms`.
Subtract those subfields to obtain the inner Schur solve time; do not add them
again to total time. The report shows this accounting explicitly. Iteration
counts sum the largest per-RHS inner count at each correction. They describe
the full Wilson operator for `none` and the half-volume Schur operator for
`even_odd`; each BiCGStab iteration needs up to two such operator applications.

This is independent batched BiCGStab, with separate coefficients and masks per
RHS. It does not share a Krylov basis. Choose batch size by time per RHS and
the measured variation, not simply by the time of a whole batch.

## Validation and remaining scope

The six-system self-test uses 2⁴ and 4⁴ volumes with RHS 1, 3, and 7, including
one zero RHS in multi-RHS tests. It checks the CPU FP64 Schur product against a
composition built from the independent direct full Wilson stencil. It also
checks transformed RHS identities, exact odd reconstruction from the known
solution, CPU and Metal FP32 products, periodic/antiperiodic boundaries,
duplicate neighbors on extent two, and rejection of odd periodic extents.
All six solver methods then undergo full FP64 residual and known-solution checks.
A separate free-field periodic solve also passes all six paths. The preexisting
nine-system unpreconditioned Wilson self-test passed as a regression after the
shared-code changes. These logs are saved in `results/evenodd-free-periodic.log`
and `results/wilson-eo-regression.log`.

The code remains a bare positive-mass Wilson prototype on synthetic SU(3) links.
It has no clover term, multigrid, domain decomposition, tuned production gauge
ensemble, or support for negative critical bare masses. No Chroma/QDP binary
is used as a numerical oracle. These results test the specified Schur formulation
and complete mixed-precision solve, not production multigrid performance.

Implementation: `evenodd.hpp` (packing, Schur operators and correction adapter),
`evenodd.metal` (combination kernel), and `evenodd.mm` (independent checks and
comparison driver). Raw timings and full-system convergence histories are in
`results/evenodd.csv` and its `.trace.csv` companion.
