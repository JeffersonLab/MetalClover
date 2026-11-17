# MetalClover API guide

This guide describes the installed public interface in
[`metalclover/metalclover.h`](../metal-bsr-prototype/include/metalclover/metalclover.h)
and its implementation in
[`metalclover_api.mm`](../metal-bsr-prototype/metalclover_api.mm).
See [BUILD.md](BUILD.md) for compilation and linking. The current header defines
`METALCLOVER_ABI_VERSION == 1` and `METALCLOVER_BATCH_API_VERSION == 1`.
Check the batch macro when building clients that require the array extension.

## Lifecycle and solver choice

1. Call `metalclover_default_options(&options)`; then set geometry and action
   parameters to match your supplied gauge state.
2. Call `metalclover_create(&options, links, &context)`. Creation copies the
   supplied gauge and builds the operator and, for MG, one hierarchy. The input
   gauge buffer may be released after creation returns.
3. Optionally set the resident RHS cap with `metalclover_set_rhs_batch`.
4. Reuse the context for scalar or array solves with the same gauge and action.
5. Call `metalclover_destroy(context)` when finished.

`METALCLOVER_BICGSTAB` selects mixed-precision odd-Schur BiCGStab;
`METALCLOVER_MG` selects right-preconditioned flexible GMRES with a three-level
MG hierarchy. The public MG implementation uses near-null/SVD fine setup and
eight retained vectors per chirality at both transfers. The standalone drivers
expose other experimental setups, including filtered fine bases; those are not
additional public API option values.

Contexts own gauge/operator/hierarchy storage. Sources, solutions, pointer
arrays and result arrays belong to the caller. Calls on one context must be
serialized, including queries and helpers; concurrent use and distributed MPI
decomposition are unsupported. Recreate a context when changing gauge, action,
mass or hierarchy parameters. Changing the RHS cap reuses its hierarchy.

In C, the declarations are global. In C++, include the header and qualify names
with `MetalCloverEnv::`, or use `using namespace MetalCloverEnv;`. The exported
symbols have C linkage, allowing GCC-built clients to call the Apple Clang dylib.

## Exact buffer layout

All public fields use interleaved FP64 real/imaginary pairs (two `double`s per
complex value). A `double` array is a portable way to supply this representation;
the API takes `void *` and does not require a particular C++ complex type.

For shape `[Lx,Ly,Lz,Lt]`, let `V=Lx*Ly*Lz*Lt` and
`site=x+Lx*(y+Ly*(z+Lz*t))`. Directions `mu=0,1,2,3` correspond to x,y,z,t.

| Buffer | Complex count | Complex index |
| --- | ---: | --- |
| Full spinor | `12*V` | `(site*3+color)*4+spin` |
| Odd spinor | `12*V/2` | `(odd_index*3+color)*4+spin` |
| Forward gauge links | `36*V` | `((site*4+mu)*3+row)*3+col` |

Real and imaginary parts of complex index `i` occupy doubles `2*i` and `2*i+1`.
Colors and matrix row/column indices range from 0 to 2; spins range from 0 to 3.
Odd fields enumerate sites with `(x+y+z+t)%2 == 1` in the same lexicographic
order. This is the library's color-then-spin layout; do not pass another
framework's raw spinor storage without checking/repacking its order and basis.

The gauge argument contains already stout-smeared links in the action's gauge
state. The library does **not** smear links. Supply the corresponding mass,
anisotropy and clover coefficients, without baking those weights into links.
Set `apply_temporal_antiperiodic=1` for un-boundary-conditioned stout links
when temporal antiperiodicity is needed; set it to 0 when links already contain
the boundary condition (or when using unconditioned periodic links). Applying
the temporal sign twice changes the operator.

Each pointer passed to a solve denotes one complete scalar-RHS field. Array
calls take arrays of field pointers; callers do not pack the GPU's internal
component-major tile format. Allocate the counts above times `2*sizeof(double)`
bytes for each field. All outputs must be distinct and disjoint from every
source field in an array call. Initialize every solution: `x` is an input guess
as well as an output. Use zeros to start from zero. Keep buffers and callback
state alive until the synchronous call returns.

## Options

Always obtain defaults from the library, rather than zero-initializing an
options struct alone. `shape` and `mass` default to zero; set them explicitly.

| Member | Default | Meaning / constraint |
| --- | --- | --- |
| `abi_version` | ABI macro | Must match the installed ABI. |
| `shape[4]` | all 0 | Positive even lattice extents, x,y,z,t. |
| `mass` | 0 | Bare mass for the supplied action; finite. |
| `xi0`, `nu` | 1, 1 | Action anisotropy parameters; positive and finite. |
| `clover_r`, `clover_t` | 1, 1 | Spatial/temporal clover coefficients; finite. |
| `apply_temporal_antiperiodic` | 1 | Apply temporal wrap sign to supplied links. |
| `method` | `METALCLOVER_MG` | MG or BiCGStab enum. |
| `workers` | 8 | CPU worker count, 1–64. |
| `restart` | 8 | Fine MG FGMRES restart, 2–12. |
| `blocks[4]` | 4,4,4,4 | Fine aggregate extents. Each must divide its lattice extent and leave an even coarse extent. |
| `null_candidates` | 24 | Fine setup candidates, 8–64 (at least retained basis size). |
| `null_iterations` | 128 | Positive setup FGMRES iteration cap per candidate. |
| `null_restart` | 8 | Setup FGMRES restart, 2–24. |
| `null_warm` | 32 | Nonnegative homogeneous warm-sweep count. |
| `null_batch` | 4 | Incremental candidate append width, 1–32; 0 keeps all candidates for a final SVD. |
| `null_tolerance` | 0.001 | Setup target, strictly between 0 and 1. Capped setup is allowed. |
| `second_filter` | 64 | Positive second-transfer filtering sweep count. |
| `bottom_iterations` | 8 | Positive bottom coarse-solver cap. |

MG-specific settings do not configure BiCGStab. The second aggregation uses
`2,2,2,2` coarse sites, hence the even first-coarse-extent requirement. For
example, a `12,12,12,96` lattice can use `3,3,3,8` fine blocks. Default `4^4`
blocks are invalid for that shape because the spatial first-coarse extents
would be odd. Creation may fail on unsuitable local clover pivots, rank-deficient
setup, unavailable Metal hardware or insufficient resources as well as invalid
parameters; these defaults are not a promise of convergence for every action.

## Full and reduced solves

```c
int metalclover_solve_full(context, b64, x64, tolerance, max_iterations,
                          apply, user, info);
int metalclover_solve_odd(context, b64, x64, tolerance, max_iterations,
                         apply, user, info);
```

These are signature summaries; use the header for typed declarations.
`tolerance` must be strictly between 0 and 1 and `max_iterations` must be
positive. It is the **total Krylov budget per RHS across FP64 defect corrections**,
not a budget reset at each refinement. The implementation also limits outer
refinement corrections to 16. Both source and solution are FP64.

The full solve accepts a full source and guess and handles parity reduction
and reconstruction. Success requires the full CPU FP64 residual norm ratio
to meet tolerance. For zero sources the norm is absolute.

The odd solve instead solves `S*x_o=b'_o`, where
`S=A_oo-H_oe*A_ee^-1*H_eo`. Its source must **already be reduced**. Its acceptance
residual is on that odd system, not the reconstructed full equation.
The following CPU FP64 helpers are provided:

| Function | Action |
| --- | --- |
| `metalclover_reduce_rhs(ctx, full_b, odd_b)` | `b'_o=b_o-H_oe*A_ee^-1*b_e` |
| `metalclover_reconstruct(ctx, full_b, odd_x, full_x)` | Copy odd solution and compute `x_e=A_ee^-1*(b_e-H_eo*x_o)` |
| `metalclover_apply_full(ctx, x, y)` | Compute `y=M*x` in CPU FP64 |
| `metalclover_apply_odd(ctx, x, y)` | Compute `y=S*x` in CPU FP64 |

Use distinct helper input/output buffers. For an odd workflow, reduce the
original full source, solve the reduced field, reconstruct with that original
source, and independently evaluate the full residual before accepting it.

## Multi-RHS arrays and actual GPU batching

```c
int metalclover_solve_full_array(ctx, rhs_count, b_ptrs, x_ptrs,
                                tolerance, max_iterations,
                                apply, user, info_array, failed_rhs);
int metalclover_solve_odd_array(ctx, rhs_count, b_ptrs, x_ptrs,
                               tolerance, max_iterations,
                               apply, user, info_array, failed_rhs);
```

For nonempty calls, supply `rhs_count` source pointers, output pointers and
`metalclover_solve_info` elements. `info_array[i]` corresponds to column `i`.
Every column has its own guess, source norm, Krylov history, convergence mask
and total iteration budget. Operator data and the hierarchy are shared. These
are independent solves in genuine GPU tiles, not a coupled block-Krylov method.
The library normalizes residuals before FP32 conversion and rescales corrections
in FP64 to support heterogeneous source magnitudes.

`metalclover_set_rhs_batch(ctx, cap)` accepts 1–64; the default is 8 and 1
selects scalar width. The cap is independent of `rhs_count`: a larger array is
split into sequential tiles, with a smaller tail. A memory-budget estimate
may lower the actual resident width. Changing widths can allocate/cache new
workspaces and pipelines without rebuilding the hierarchy. Bound caller-side
groups as well: a GPU cap does not eliminate memory occupied by the caller's
full FP64 spinors.

After a call, use `metalclover_get_batch_info(ctx, &batch)` to inspect actual
`maximum_resident_rhs` and `tiles`; do not infer width from array length alone.
An empty array succeeds with a valid context and need not supply fields/info.
`failed_rhs` is optional; when provided it is `SIZE_MAX` on success. On failure,
it identifies the first unaccepted column of a failed tile, or the precise
column for validation/callback failures. Read the status first. Preserve
columns whose `info[i].converged` is true; other outputs may contain partial
iterates and are not accepted solutions. Input-validation failures occur before
GPU execution. Do not blindly retry every source or assume failure rolled back
the entire array.

For distillation, one color eigenvector normally supplies **four source-spin
RHS**. Thus 24 color vectors correspond to 96 inversions, not 24. Source building
and sink projection are caller/Chroma operations, not C ABI array packing.

## FP64 callback and errors

The optional callback type is:

```c
typedef int (*metalclover_apply64)(void *user, const void *input,
                                 void *output, size_t n_complex);
```

It must overwrite `output` with the **PLUS-sign** action `M*input` for full
solves, or `S*input` for odd solves, in the same FP64 layout and subset.
It receives one scalar field at a time even in an array solve. Return 0 on
success; any other return fails the solve. The same `user` is forwarded on
each call. Match exactly the gauge, basis, action and boundary conventions
used at creation. Do not recursively call context operations from the callback.
Pass `nullptr, nullptr` (C: `NULL, NULL`) to use the library's CPU FP64 operator.
An independent callback, for example original Chroma, lets the outer refinement
check against the application operator.

Integer-returning API functions return 0 on success and nonzero on failure.
Read `metalclover_last_error()` immediately on the calling thread; copy it if
needed because a later successful protected call clears the thread-local error.
`metalclover_create` clears its output context on entry, including failures.
Destroy successful contexts exactly once. Required output arguments such as
solve info must be non-null. `default_options` and `destroy` return void;
`setup_seconds` returns a time rather than a status.

## Result fields

`metalclover_solve_info` reports `converged`, `relative_residual`,
`outer_iterations`, `krylov_iterations`, `fine_schur_applications`,
`seconds`, `residual_seconds`, and nominal operation counts/rates. Work fields
are `nominal_operator_flops`, `fine_flops`, `coarse_flops`, `transfer_flops`,
and `local_vector_flops`. Array column times are amortized shares of the total
array-call time, not independently measured column latencies. Dispatched work,
including inactive lanes, is apportioned over corrected columns.

`metalclover_batch_info` reports `rhs_count`, `configured_batch_rhs`,
`maximum_resident_rhs`, `tiles`, `seconds`, `workspace_seconds`,
`residual_seconds`, `nominal_operator_flops`, `fine_schur_dispatches`, and
`nominal_operator_gflops`. Use aggregate batch time/work/rate for throughput.
Fine Schur dispatches are not per-column Krylov iterations. First-use batch
workspace creation is included in call time and also reported separately.

`metalclover_setup_seconds(ctx)` reports context setup separately. Nominal
GFLOP/s counts source-level dispatched operator work, including padded coarse
work; it excludes Krylov-vector and CPU FP64 work from the numerator while
including their time in the denominator. It is not a hardware-counter rate.
See [the detailed accounting](../metal-bsr-prototype/OPERATOR_WORK.md).

## C++ multi-RHS example

This function accepts caller-prepared stout links and full spinor sources,
allocates zero guesses, submits one array call, and releases the context on
all exit paths. It can be compiled as a translation unit; supply your own
input reader and `main`. Keep/reuse a context outside such a wrapper for repeated
calls so hierarchy setup is amortized. The example throws on API failure
instead of accepting partial results.

```cpp
#include <metalclover/metalclover.h>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

using namespace MetalCloverEnv;

std::vector<std::vector<double>> solve_sources(
    const metalclover_options &options,
    const std::vector<double> &stout_links,
    const std::vector<std::vector<double>> &sources) {
  size_t volume = 1;
  for (auto extent : options.shape) volume *= extent;
  const size_t spinor_doubles = 24 * volume;
  if (stout_links.size() != 72 * volume)
    throw std::runtime_error("wrong gauge size");
  for (const auto &b : sources)
    if (b.size() != spinor_doubles)
      throw std::runtime_error("wrong source size");

  auto check = [](int status) {
    if (status) throw std::runtime_error(metalclover_last_error());
  };
  metalclover_context *raw = nullptr;
  check(metalclover_create(&options, stout_links.data(), &raw));
  std::unique_ptr<metalclover_context, decltype(&metalclover_destroy)>
      ctx(raw, metalclover_destroy);
  check(metalclover_set_rhs_batch(ctx.get(), 8));

  const size_t count = sources.size();
  std::vector<std::vector<double>> solutions(
      count, std::vector<double>(spinor_doubles, 0.0));
  std::vector<const void *> b(count);
  std::vector<void *> x(count);
  std::vector<metalclover_solve_info> info(count);
  for (size_t i = 0; i < count; ++i) {
    b[i] = sources[i].data();
    x[i] = solutions[i].data();
  }
  size_t failed = SIZE_MAX;
  int status = metalclover_solve_full_array(
      ctx.get(), count, b.data(), x.data(), 1e-12, 20000,
      nullptr, nullptr, info.data(), &failed);
  if (status) {
    // Copy the error before any other API query clears it.
    std::string error = metalclover_last_error();
    throw std::runtime_error("array failed at column " +
                             std::to_string(failed) + ": " + error);
  }
  for (const auto &result : info)
    if (!result.converged)
      throw std::runtime_error("unaccepted result");
  metalclover_batch_info batch{};
  check(metalclover_get_batch_info(ctx.get(), &batch));
  std::printf("%zu RHS, max resident %zu, %zu tiles, %.6f s\n",
              batch.rhs_count, batch.maximum_resident_rhs,
              batch.tiles, batch.seconds);
  return solutions;
}
```

Before calling this function, initialize options with `default_options`, then
set the actual action. For an isotropic 8^4 example, set all four shape entries
to 8 and set `mass=0.1`; the remaining action defaults are isotropic. This is
only an example configuration, not parameters for an arbitrary physical field.
For large arrays, call a reusable context with bounded groups rather than
allocating all sources/solutions as this small example does.
