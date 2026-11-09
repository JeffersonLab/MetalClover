#!/usr/bin/env python3
"""Validate raw sparse-MG measurements and regenerate the 16^4 report."""
import argparse
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

ap = argparse.ArgumentParser()
ap.add_argument("csv", nargs="?", type=Path, default=Path("results/sparse-mg16.csv"))
ap.add_argument("--output", type=Path, default=Path("SPARSE_MULTIGRID_RESULTS.md"))
args = ap.parse_args()
with args.csv.open() as f:
    rows = list(csv.DictReader(f))
groups = defaultdict(list)
for r in rows:
    assert all(math.isfinite(float(r[k])) for k in ("total_ms", "max_true_residual", "solution_error", "worst_coarse_true_residual"))
    assert float(r["max_true_residual"]) <= float(r["tolerance"]) == 1e-10
    assert float(r["solution_error"]) <= 1e-7
    assert int(r["extent"]) == 16 and float(r["csw"]) == 1
    assert int(r["tile_rhs"]) == 96 and int(r["rhs"]) == 96*int(r["tiles"])
    assert int(r["coarse_dim"]) == 6144 and int(r["coarse_nnz_blocks"]) == 8448
    assert int(r["coarse_limit"]) == 4 and int(r["vectors"]) == 12
    groups[float(r["mass"]), int(r["rhs"]), r["backend"]].append(r)
cases = sorted({k[:2] for k in groups})
assert cases == [(.1,96),(.1,192),(.5,96),(.5,192)]
for k in cases:
    for b in ("cpu_sparse_mg", "metal_sparse_mg"):
        assert len(groups[k + (b,)]) == 3, (k,b)


def values(k, b, field):
    return [float(r[field]) for r in groups[k + (b,)]]


def med(k, b, field="total_ms"):
    return st.median(values(k,b,field))


def span(k, b, field):
    v = sorted({int(r[field]) for r in groups[k+(b,)]})
    return str(v[0]) if len(v) == 1 else f"{v[0]}–{v[-1]}"


lines = ["# Sparse coarse multigrid: 16⁴ CPU versus Metal", "",
         "Local measurements on Apple M1 Max (24 GPU cores, 32 GB, macOS "
         "26.7.1). The standalone clover-Wilson operator keeps c_sw=1, QDP "
         "gamma conventions, gauge strength 0.25, positive bare masses "
         "0.1 and 0.5, and antiperiodic time. The comparison uses the same "
         "sparse two-level solver on CPU and Metal, with complete CPU FP64 "
         "residual checks on the reconstructed full equation for every RHS.", "",
         "## Complete solve time", "",
         "Medians of three measured logical batches, from zero. **96 RHS "
         "are resident simultaneously. 192 RHS run as two sequential "
         "96-RHS tiles**, sharing the hierarchy/operator/inner workspace. CPU and "
         "Metal workspaces are allocated sequentially to fit the 32 GB Mac. "
         "One complete resident-tile warmup precedes each backend's timings "
         "for each mass. Backend order reverses between masses. CPU uses "
         "eight GCD workers and single-threaded Accelerate for setup. "
         "The computer had ordinary local activity, without thermal or "
         "scheduler isolation.", "",
         "| Mass | Logical RHS | Resident RHS | CPU total s | Metal total s | CPU / Metal | CPU RHS/s | Metal RHS/s |",
         "|---:|---:|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    c,g = (med(k,b) for b in ("cpu_sparse_mg","metal_sparse_mg"))
    lines.append(f"| {k[0]:g} | {k[1]} | 96 | {c/1000:.3f} | {g/1000:.3f} | {c/g:.2f}× | {1000*k[1]/c:.2f} | {1000*k[1]/g:.2f} |")
lines += ["", "Solve totals include resets, FP32 RHS reduction and even-site "
          "reconstruction, fine and coarse sparse products, smoothing, "
          "FGMRES/BiCGStab vector operations and reductions, CPU small QR, "
          "GPU submissions/waits, FP64 accumulation, and full FP64 residual "
          "checks. Synthetic manufactured RHS generation is outside solve "
          "time and has its own raw column; supplied RHS would not require "
          "that operator application. Hierarchy, problem, and inner workspace "
          "setup are reported separately. FP64 outer buffers are allocated "
          "once per logical case outside the solve timer; their allocation "
          "time was not collected. Reset and data movement inside solve are "
          "timed. No "
          "192-RHS result represents 192 simultaneously resident Krylov vectors.", "",
          "| Mass | RHS | CPU min–max s | Metal min–max s | Metal CPU-FP64 residual s | Metal FP32 inner s |",
          "|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    c,g = (values(k,b,"total_ms") for b in ("cpu_sparse_mg","metal_sparse_mg"))
    lines.append(f"| {k[0]:g} | {k[1]} | {min(c)/1000:.3f}–{max(c)/1000:.3f} | {min(g)/1000:.3f}–{max(g)/1000:.3f} | {med(k,'metal_sparse_mg','residual_ms')/1000:.3f} | {med(k,'metal_sparse_mg','inner_ms')/1000:.3f} |")

lines += ["", "## Convergence and work", "",
          "The fine solver is FGMRES(8) with an inner target of 0.002. The "
          "basis has twelve vectors per chirality, 4⁴ aggregates, eight "
          "FP64 filtering steps at weight 0.65, and one pre/post clover-block "
          "Jacobi step at weight 1.0. Fine S counts include smoothing, "
          "Arnoldi applications, and actual FP32 residual refreshes. Counts "
          "are summed over resident tiles for the logical batch; a call "
          "always applies to 96 resident RHS.", "",
          "| Mass | RHS | Backend | Outer corrections/tile | Fine Krylov steps/tile | Total fine S calls | Total coarse SpMV dispatches | Total coarse iterations | Worst full FP64 residual |",
          "|---:|---:|---|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    for b in ("cpu_sparse_mg","metal_sparse_mg"):
        lo = min(values(k,b,"krylov_min")); hi = max(values(k,b,"krylov_max"))
        outerlo = min(values(k,b,"outer_min")); outerhi = max(values(k,b,"outer_max"))
        kry = f"{lo:g}" if lo == hi else f"{lo:g}–{hi:g}"
        outer = f"{outerlo:g}" if outerlo == outerhi else f"{outerlo:g}–{outerhi:g}"
        lines.append(f"| {k[0]:g} | {k[1]} | {b} | {outer} | {kry} | {span(k,b,'fine_schur_calls')} | {span(k,b,'coarse_spmv_calls')} | {span(k,b,'coarse_iterations')} | {max(values(k,b,'max_true_residual')):.3g} |")
lines += ["", "The coarse solver is left block-Jacobi-preconditioned "
          "complex BiCGStab, starting from zero, with target relative "
          "preconditioned residual 0.1 and a four-iteration cap. The returned "
          "coarse correction may be inexact; FGMRES stores Z to accommodate "
          "that variable preconditioner. The unpreconditioned FP32 coarse "
          "residual is recomputed after every correction. It is a diagnostic "
          "and can exceed 0.1. CPU stops submitting sparse products after "
          "convergence; Metal encodes the cap with per-RHS masking. GPU "
          "SpMV counts include masked/padded dispatches and both backends "
          "include one true-coarse-residual product per correction. This "
          "difference is included in the runtime comparison.", "",
          "| Mass | RHS | Worst CPU coarse true residual | Worst Metal coarse true residual |",
          "|---:|---:|---:|---:|"]
for k in cases:
    lines.append(f"| {k[0]:g} | {k[1]} | {max(values(k,'cpu_sparse_mg','worst_coarse_true_residual')):.4f} | {max(values(k,'metal_sparse_mg','worst_coarse_true_residual')):.4f} |")

with Path(str(args.csv)+".trace.csv").open() as f:
    traces = list(csv.DictReader(f))
for r in traces:
    assert math.isfinite(float(r["max_true_residual"]))
history = defaultdict(dict)
for r in traces:
    if float(r["mass"]) == .1 and r["rhs"] == "192" and r["trial"] == "0":
        history[int(r["offset"]),int(r["outer_step"])][r["backend"]] = float(r["max_true_residual"])
lines += ["", "Example full FP64 convergence history, mass 0.1, 192 RHS, trial 0:", "",
          "| Global RHS offset | Outer correction | CPU full residual | Metal full residual |",
          "|---:|---:|---:|---:|"]
for (offset,step), v in sorted(history.items()):
    lines.append(f"| {offset} | {step} | {v['cpu_sparse_mg']:.3g} | {v['metal_sparse_mg']:.3g} |")

lines += ["", "## Sparse storage, setup, and memory", "",
          "At 16⁴ the retained odd space has 393,216 complex components. "
          "256 aggregates × 24 components give 6,144 coarse unknowns. The "
          "exact Galerkin operator has 8,448 dense 24×24 blocks (33 per "
          "aggregate), stored as block CSR. Its FP32 coefficients occupy "
          "**37.125 MiB**, compared with 288 MiB for a full FP32 coarse "
          "matrix. An FP64 sparse copy, sparse indices, basis, and local "
          "diagonal inverses are additional. Only local 24×24 blocks are "
          "inverted; no global dense inverse, factorization, or fine-by-all-"
          "coarse prolongator array is formed. The latter would take 36 GiB "
          "per FP64 array at this geometry, before all other workspaces.", "",
          "| Mass | RHS | Hierarchy setup s | CPU inner workspace setup s | Metal inner workspace setup s | Problem setup s | Process peak RSS GiB |",
          "|---:|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    peak = max(values(k,'cpu_sparse_mg','peak_rss_mib') + values(k,'metal_sparse_mg','peak_rss_mib'))/1024
    lines.append(f"| {k[0]:g} | {k[1]} | {med(k,'cpu_sparse_mg','hierarchy_setup_ms')/1000:.3f} | {med(k,'cpu_sparse_mg','workspace_setup_ms')/1000:.3f} | {med(k,'metal_sparse_mg','workspace_setup_ms')/1000:.3f} | {med(k,'cpu_sparse_mg','problem_setup_ms')/1000:.3f} | {peak:.2f} |")
lines += ["", "Hierarchy setup includes FP64 test-vector filtering, local "
          "QR, sparse assembly, local block inversion, and temporary setup "
          "operator construction. One mass-specific hierarchy is reused "
          "across both processors, 96/192 RHS, and all trials. Process RSS "
          "is the cumulative OS high-water mark within this run; it is not "
          "a per-method isolated measurement or total system memory usage. "
          "Unrestricted simultaneous 192-RHS fine Krylov stacks would impose "
          "a much larger memory footprint. The user-visible 192-RHS batch "
          "therefore uses two bounded-memory tiles. With many RHS, hierarchy "
          "setup is paid once per gauge/mass/geometry, not once per tile. "
          "The setup table is not a complete cold-start total because "
          "FP64 outer-buffer allocation was not collected.", "",
          "## Validation and reproduction", "",
          f"All {len(rows)} measured logical batches passed the full FP64 "
          "criterion for every RHS. They contain "
          f"{sum(int(r['tiles']) for r in rows)} resident-tile solves and "
          f"{sum(int(r['rhs']) for r in rows)} RHS executions across repeated "
          "trials. Worst final true relative residual was "
          f"{max(float(r['max_true_residual']) for r in rows):.3g}; maximum "
          "known-solution relative error was "
          f"{max(float(r['solution_error']) for r in rows):.3g}. The trace "
          "records each tile's global RHS offset and every "
          "outer true residual. A known solution is used only for validation "
          "and RHS generation; hierarchy setup is independent of it.", "",
          "Before timings, the 40-solve clover regression, 24-solve dense "
          "multigrid regression, and eight-solve sparse self-test passed. "
          "Sparse self-tests compare with independent dense column probing, "
          "coarse gamma5 Hermiticity, CPU/Metal sparse products, and strict "
          "coarse solves, including zero RHS locking. Actual 8⁴ and 16⁴ "
          "hierarchies passed FP64 Galerkin/coarse-solver checks at both "
          "masses. At 16⁴ the Galerkin random-vector identity agrees at "
          "about 2e-15 relative error.", "",
          "```sh", "make build/metal_sparse_multigrid build/metal_multigrid build/metal_clover",
          "python3 run_sparse_study.py", "python3 sparse_summary.py", "```", "",
          "The main timing command alone is:", "", "```sh",
          "VECLIB_MAXIMUM_THREADS=1 ./build/metal_sparse_multigrid \\",
          "  --sizes 16 --rhs 96,192 --tile-rhs 96 --mass .1,.5 --csw 1 \\",
          "  --vectors 12 --aggregate 4 --restart 8 --coarse-iters 4 \\",
          "  --coarse-tol .1 --samples 3 --threads 8 --csv results/sparse-mg16.csv",
          "```", "",
          "This establishes CPU/Metal convergence and timing for the "
          "sparse iterative coarse solver at the larger volume. The gauge "
          "field remains synthetic and the masses positive; production "
          "ensembles, critical-mass behavior, distributed execution, and "
          "additional multigrid levels need separate validation. "
          "See [SPARSE_MULTIGRID.md](SPARSE_MULTIGRID.md) for the assembly "
          "equations and implementation details.", "",
          f"Raw data: [{args.csv.as_posix()}]({args.csv.as_posix()}), "
          f"[full FP64 traces]({args.csv.as_posix()}.trace.csv), "
          "[benchmark log](results/sparse-mg16.log), "
          "[regression/study log](results/sparse-study.log).", ""]
args.output.write_text("\n".join(lines))
print("Wrote", args.output, ":", len(rows), "logical batches")
