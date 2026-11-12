# Operator-work accounting and solver interfaces

All C++ library types and routines reside in the unique `MetalCloverEnv`
namespace, including physics, Krylov, transfer and storage helpers. Standalone
programs import that namespace explicitly. The prefixed `metalclover_*` entry points have C linkage for the
cross-compiler ABI; their C++ declarations also reside in `MetalCloverEnv`.
C callers use the same unmangled prefixed symbols.

The solver returns a `metalclover_solve_info` through the public C ABI in
`include/metalclover/metalclover.h`. The caller chooses how to print it. Both
inverters report iterations, full or Schur relative residual, elapsed seconds,
FP64 residual-check seconds, fine Schur applications and nominal operator work.
The standalone low-memory driver writes the same work totals and rate to CSV.
Hierarchy construction is separate from each reused-hierarchy solve.

## What the rate means

`nominal_operator_gflops = nominal_operator_flops / (seconds * 1e9)`.
This is a source-level operator-work model, **not hardware-counter FLOP/s**.
A real multiply or add counts as one operation; FMA counts as two; a complex
multiply/add pair counts as eight. Division counts as one. Signs, conjugation,
integer addressing, loads, stores and copies count as zero. Explicit
component multiplications by constants (including the phase sign scale)
are counted before compiler simplification; unary negation is not.

The numerator includes dispatched fine hopping and clover kernels, coarse
sparse products and diagonal applications, restriction/prolongation and the
local vector arithmetic adjoining these operators. Coarse BiCGStab uses a
fixed dispatch limit with device masks; its padded products are counted even
when a converged RHS masks their arithmetic. Thus this part is nominal work,
not a measurement of executed instructions. The numerator excludes Krylov
vector updates, dot reductions, scalar QR, CPU FP64 residual evaluation and
CPU clover regeneration. It must not be described as total solver FLOP/s or
compared directly with peak hardware throughput.

The denominator is the complete solve time including refinement, solution
updates and CPU FP64 checks. In the low-memory benchmark it includes FP64
solution slab writeback and gauge/residual page management. The final
independent native-Chroma verification process is reported separately.
Chroma's adapter adds field conversion and its final native residual check to
its reported solve time. Setup/SVD work is excluded from the solve rate.

## Formulas

For each complex output of `compact_hop_rect`, each of eight directions
performs two real phase scales plus, for each of three colors, two gauge
scales, two complex multiplies and two complex additions: `8*(2+3*18)=448`.
There are 12 complex outputs per destination site: 5376 operations/site.

The current scalar-output LDL kernel repeats a complete six-component
triangular solve for each requested output: `(15+15)*8+12=252`
operations/output, or 3024/site. Its matrix multiply takes 152/output on
average, or 1824/site. Combining it with the hopping product costs a further
four operations/output. A finest odd Schur product therefore counts
`2*5376+3024+1824+48=15648` operations/odd site/RHS. This describes the current
kernels; it is deliberately not the conventional optimized Wilson Dslash
1320-FLOP convention.

A dense/coarse complex matrix multiply counts eight operations per matrix
entry/RHS, including structural entries stored inside each sparse block.
Diagonal kernels use their stored block dimension. Restriction and
prolongation each count `8 * fine_complex_dofs * vectors_per_chirality * RHS`;
prolongation also adds two operations per fine complex component for its
accumulation into the solution. Complex differences/additions count two per
component; real scale and real AXPY count two and four respectively.

The returned fine/coarse/transfer/local-vector subtotals make the scope
inspectable. `fine_flops` counts the compact fine-operator and reduction/
reconstruction kernels; `coarse_flops` groups coarse matrix kernels and MG
diagonal smoothers at every level, including the fine clover smoother. Large coarse-padding contributions should be interpreted using
the separately logged actual and padded bottom iterations.

## Full and odd interfaces

`metalclover_solve_full` accepts a full source and initial guess. Its FP32
corrections perform even/odd reduction and reconstruction on Metal; the
accumulated solution and every outer residual use FP64.

`metalclover_solve_odd` accepts an **already reduced odd source** for
`S_oo = A_oo - D_oe A_ee^-1 D_eo`. It returns only the odd solution and does not
reduce the RHS again. `metalclover_reduce_rhs` and `metalclover_reconstruct`
provide independent CPU FP64 conversion helpers. Reconstruction uses the
full source, not the reduced odd source.

Each context can serve either interface and reuses its selected inverter and
hierarchy. One RHS is resident per call. The C ABI has no C++ standard-library
objects or Objective-C types, so GCC/libstdc++ Chroma can link the Apple
Clang/libc++ backend. The installed dylib embeds the Metal shader sources and
does not depend on the working directory.

The low-memory production benchmark continues to stream large FP64 fields
through file-backed storage. The first C ABI/Chroma harness uses in-memory
FP64 source, solution and residual workspaces; its memory footprint in a
large Chroma process also includes Chroma's gauge and clover allocations.
The standalone 32-cubed memory result does not establish that the complete
32-cubed Chroma process fits the same 32 GiB laptop.

The `metalclover_solve_full_array` and `metalclover_solve_odd_array` C ABI calls
accept arrays of field pointers and return an info array in the same order.
They reuse one context/hierarchy and one GPU workspace, processing columns
sequentially. They report a failed-column index; an empty array succeeds.
Caller-owned output fields must be distinct and disjoint from source fields.
The native Chroma overrides use the corresponding vector-of-shared-pointers
signature used by `SB::doInversion` and expose per-column info plus an array
operator-work rate. The caller still owns storage for the arrays and can
stream smaller batches from a much larger RHS collection.

The dylib defaults to a macOS 15.0 deployment target. Override
`METALCLOVER_MACOS_MIN` when building for another supported macOS release;
Metal device and shader support are checked when the context is created.

Runtime shader compilation explicitly selects Metal 3.1, so a GCC-built
calling executable cannot silently select the legacy Metal 1.1 default.
