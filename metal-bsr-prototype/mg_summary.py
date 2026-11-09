#!/usr/bin/env python3
"""Report complete CPU/Metal two-level clover solve measurements."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("--output", type=Path, default=Path("MULTIGRID_RESULTS.md"))
    args = ap.parse_args()
    with args.csv.open() as f:
        data = list(csv.DictReader(f))
    if not data:
        raise SystemExit("No measurements")
    groups = defaultdict(list)
    for row in data:
        if float(row["max_true_residual"]) > float(row["tolerance"]):
            raise SystemExit("Full FP64 clover residual target missed")
        key = (float(row["mass"]), int(row["extent"]), int(row["rhs"]))
        groups[key + (row["backend"],)].append(row)
    cases = sorted({k[:3] for k in groups})
    names = ("cpu_bicgstab", "metal_bicgstab", "cpu_fgmres_smooth", "metal_fgmres_smooth", "cpu_fgmres_mg", "metal_fgmres_mg")
    for key in cases:
        for name in names:
            if len(groups[key + (name,)]) != 3:
                raise SystemExit(f"Expected three trials for {key}, {name}")

    def med(key, name, field="total_ms"):
        return statistics.median(float(r[field]) for r in groups[key + (name,)])

    def counts(key, name):
        values = []
        for field in ("outer_steps", "inner_max_sum", "fine_schur_calls", "coarse_calls"):
            v = sorted({int(r[field]) for r in groups[key + (name,)]})
            values.append(str(v[0]) if len(v) == 1 else f"{v[0]}–{v[-1]}")
        return " / ".join(values)

    lines = ["# Two-level clover multigrid on M1 Max", "",
             "Measured locally on November 4, 2026, on the M1 Max (24 GPU cores, "
             "32 GB, macOS 26.7.1, AC power). CPU runs use eight worker tasks "
             "and single-threaded Accelerate for the small coarse solve. "
             "Entries are medians of three complete solves from zero after a "
             "full warmup, with rotating method order. All methods share the "
             "same c_sw=1 clover operator and RHS within a case.", "",
             "**The two-level implementation passed its numerical checks, "
             "but was slower than the existing BiCGStab baseline in these tests.** "
             "The coarse grid reduces FGMRES iteration counts; the tables below "
             "include the extra fine Schur applications and transfer work needed "
             "by each cycle, rather than judging performance by iterations alone.", "",
             "## Complete solve times", "",
             "BiCGStab is the existing even/odd clover baseline with FP32 "
             "elimination. FGMRES uses the same odd Schur operator and FP32 "
             "elimination; its preconditioner is either two clover-block smoothing "
             "steps or the two-level cycle. Full CPU FP64 residual checks and "
             "FP64 accumulation are retained for every method. Times include "
             "preparation, reconstruction, vector operations, reductions, small "
             "FGMRES QR solves, GPU submission/waits, and FP64 checks. Setup is "
             "excluded from warm solve times and reported separately.", "",
             "| Mass | RHS | CPU BiCGStab (ms) | CPU two-level (ms) | Metal BiCGStab (ms) | Metal two-level (ms) | Metal MG gain vs BiCGStab | CPU MG / Metal MG |",
             "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        cb, cm, gb, gm = (med(key, b) for b in ("cpu_bicgstab", "cpu_fgmres_mg", "metal_bicgstab", "metal_fgmres_mg"))
        lines.append(f"| {m:g} | {n} | {cb:.2f} | {cm:.2f} | {gb:.2f} | {gm:.2f} | {gb/gm:.2f}× | {cm/gm:.2f}× |")
    lines += ["", "Gain is baseline time divided by multigrid time; below one "
              "means multigrid is slower. These comparisons use fresh timings "
              "from the same run, not historical clover results.", "",
              "## Isolating the coarse-grid effect", "",
              "| Mass | RHS | CPU smooth FGMRES (ms) | CPU two-level (ms) | Metal smooth FGMRES (ms) | Metal two-level (ms) | Metal coarse-grid gain |",
              "|---:|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        cs, cm, gs, gm = (med(key, b) for b in ("cpu_fgmres_smooth", "cpu_fgmres_mg", "metal_fgmres_smooth", "metal_fgmres_mg"))
        lines.append(f"| {m:g} | {n} | {cs:.2f} | {cm:.2f} | {gs:.2f} | {gm:.2f} | {gs/gm:.2f}× |")
    lines += ["", "The smoother-only path has two damped block-Jacobi steps "
              "and one fine Schur application per preconditioner call. The "
              "two-level path inserts a coarse correction and recomputes the "
              "residual before post-smoothing, using two fine Schur applications "
              "per call. FGMRES then needs its own Schur application. A coarse-grid "
              "iteration reduction must offset that added cost.", "",
              "## Convergence and actual operator work", "",
              "Counts are `outer FP64 corrections / summed max per-RHS Krylov "
              "iterations / actual fine Schur calls / coarse calls`. BiCGStab "
              "and FGMRES iterations represent different work. Metal BiCGStab "
              "includes command-chunk padding; FGMRES has no such padding. "
              "FGMRES recomputes its true FP32 Schur residual at every restart, "
              "including when its small QR estimate suggests convergence.", "",
              "| Mass | RHS | Method | Corrections / iterations / fine S / coarse | Worst full FP64 residual |",
              "|---:|---:|---|---:|---:|"]
    for key in cases:
        m, l, n = key
        for name in names:
            worst = max(float(r["max_true_residual"]) for r in groups[key + (name,)])
            lines.append(f"| {m:g} | {n} | {name} | {counts(key,name)} | {worst:.3g} |")
    with Path(str(args.csv)+".trace.csv").open() as f:
        traces = list(csv.DictReader(f))
    example = min(cases, key=lambda k:(k[0],-k[1],-k[2]))
    history = defaultdict(dict)
    for tr in traces:
        if (float(tr["mass"]),int(tr["extent"]),int(tr["rhs"])) == example and tr["trial"] == "0" and tr["backend"].startswith("metal"):
            history[int(tr["outer_step"])][tr["backend"]] = float(tr["max_true_residual"])
    lines += ["", f"Full-system residual history for mass {example[0]:g}, {example[2]} RHS:", "",
              "| Correction | Metal BiCGStab | Metal smoother FGMRES | Metal two-level FGMRES |", "|---:|---:|---:|---:|"]
    for step, vals in sorted(history.items()):
        lines.append(f"| {step} | {vals.get('metal_bicgstab',float('nan')):.3g} | {vals.get('metal_fgmres_smooth',float('nan')):.3g} | {vals.get('metal_fgmres_mg',float('nan')):.3g} |")
    lines += ["", "## Hierarchy and amortization", "",
              "The 8⁴ fine lattice has 24,576 complex odd-site components. "
              "4⁴ aggregation leaves 16 aggregates with 12 coarse components "
              "each (six per chirality), giving a 192×192 coarse operator. "
              "Six seeded full-spin test vectors with constant and random parts "
              "are filtered by eight homogeneous block-Jacobi steps, split by "
              "chirality, and orthonormalized twice locally. Setup uses no "
              "benchmark RHS or known solution. P and restriction P† preserve "
              "gamma5. The FP64 Galerkin matrix is `A_c=P† S_o P`; its direct "
              "inverse is computed on the CPU once and rounded for FP32 CPU/GPU "
              "applications. One pre- and one post-smoothing step use omega=0.65. "
              "The same mass-specific hierarchy is reused at 96 and 192 RHS.", "",
              "| Mass | RHS | Hierarchy setup (ms) | CPU MG warm ms/RHS | Metal MG warm ms/RHS | Metal MG RHS/s | Hierarchy break-even batches vs Metal BiCGStab |",
              "|---:|---:|---:|---:|---:|---:|---|"]
    for key in cases:
        m, l, n = key
        setup = med(key,"metal_fgmres_mg","hierarchy_setup_ms")
        cb, cm, gb, gm = (med(key,b) for b in ("cpu_bicgstab","cpu_fgmres_mg","metal_bicgstab","metal_fgmres_mg"))
        be = f"{setup/(gb-gm):.1f}" if gb>gm else "No warm-time gain"
        lines.append(f"| {m:g} | {n} | {setup:.2f} | {cm/n:.3f} | {gm/n:.3f} | {1000*n/gm:.1f} | {be} |")
    lines += ["", "Hierarchy setup includes filtering, local QR, probing P's "
              "columns, and dense coarse inversion. Problem generation and "
              "per-method workspace/pipeline setup are separate raw columns. "
              "BiCGStab and smoother-only solves do not need hierarchy setup; "
              "the field is repeated in their CSV rows only for case provenance. "
              "If the warm multigrid solve is slower, an unlimited supply of RHS "
              "cannot amortize away that per-solve overhead.", "",
              "## Timing spread and FP64 cost", "",
              "| Mass | RHS | CPU MG min–max (ms) | Metal MG min–max (ms) | Metal MG inner (ms) | Metal MG CPU FP64 checks (ms) |",
              "|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        ct = [float(r["total_ms"]) for r in groups[key + ("cpu_fgmres_mg",)]]
        gt = [float(r["total_ms"]) for r in groups[key + ("metal_fgmres_mg",)]]
        lines.append(f"| {m:g} | {n} | {min(ct):.2f}–{max(ct):.2f} | {min(gt):.2f}–{max(gt):.2f} | {med(key,'metal_fgmres_mg','inner_ms'):.2f} | {med(key,'metal_fgmres_mg','residual_ms'):.2f} |")
    lines += ["", "Inner time includes FP32 RHS reduction/reconstruction. "
              "The Metal FGMRES vector, transfer, coarse-solve and fine-operator "
              "work stays on the GPU; only small Hessenberg coefficients/norms "
              "are read per Arnoldi step for CPU FP64 Givens rotations and "
              "triangular solves. Those synchronizations and calculations are "
              "inside the timing. The run used ordinary local system activity, "
              "without thermal or scheduler isolation.", "",
              "## Validation and scope", "",
              f"All {len(data)} measured solves met the full CPU FP64 true "
              "relative residual target of 1e-10 for every RHS, after even-site "
              "reconstruction. Maximum known-solution relative error was "
              f"{max(float(r['solution_error']) for r in data):.3g}. "
              "Small self-tests verify P†P, coarse gamma5 Hermiticity, the "
              "Galerkin identity, exact coarse correction on range(SP), "
              "CPU/Metal restriction and prolongation, V-cycle agreement, "
              "zero RHS locking, and all six complete solver paths. The actual "
              "8⁴/4⁴-aggregate hierarchies also passed those operator identities "
              "at both benchmark masses, and the original 40-solve clover "
              "regression passed after the shared wrapper change.", "",
              "This is a first two-level implementation on small positive-mass "
              "synthetic gauge fields. It has six local vectors per chirality "
              "and fixed smoothing/filtering parameters, not a tuned production "
              "adaptive near-null solve. Dense coarse inversion is intentionally "
              "limited to 512 coarse unknowns. Different aggregates, basis "
              "sizes, smoother/coarse accuracy and setup quality need separate "
              "measurements. The result does not establish production multigrid "
              "scaling or behavior near critical mass.", "",
              "See [MULTIGRID.md](MULTIGRID.md) for equations, algorithms and reproduction.", "",
              f"Raw data: [{args.csv.as_posix()}]({args.csv.as_posix()}), "
              f"[FP64 convergence traces]({args.csv.as_posix()}.trace.csv).", ""]
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output}: {len(cases)} cases / {len(data)} solves")


if __name__ == "__main__":
    main()
