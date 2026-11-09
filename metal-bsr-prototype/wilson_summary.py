#!/usr/bin/env python3
"""Summarize the Wilson full-solve benchmark without third-party packages."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("--output", type=Path, default=Path("WILSON_RESULTS.md"))
    args = ap.parse_args()
    with args.csv.open() as f:
        data = list(csv.DictReader(f))
    if not data:
        raise SystemExit("No results")
    groups = defaultdict(list)
    for r in data:
        assert float(r["max_true_residual"]) <= float(r["tolerance"])
        groups[float(r["mass"]), int(r["extent"]), int(r["rhs"]), r["backend"]].append(r)
    cases = sorted({(m, l, n) for m, l, n, b in groups})
    backends = ("cpu_fp64", "cpu_mixed", "metal_mixed")
    for key in cases:
        for b in backends:
            if len(groups[key + (b,)]) != 3:
                raise SystemExit(f"Expected three trials for {key} / {b}")

    def med(key, b, field="total_ms"):
        return statistics.median(float(r[field]) for r in groups[key + (b,)])

    def integer_range(key, b, field):
        vals = sorted({int(r[field]) for r in groups[key + (b,)]})
        return str(vals[0]) if len(vals) == 1 else f"{vals[0]}–{vals[-1]}"

    lines = ["# Wilson / BiCGStab results on M1 Max", "",
             "Measured locally on November 3, 2026: M1 Max, 24 GPU cores, "
             "32 GB memory, macOS 26.7.1, AC power. CPU runs request eight GCD "
             "workers and use ARM NEON. Each entry is the median of three "
             "complete solves from a zero initial guess, after a full warmup. "
             "Backend order rotates between trials.", "",
             "This is the **unpreconditioned Wilson operator**, not the earlier "
             "SPD Laplacian: `(m+4)I - 0.5 Dslash`, with spin projectors, "
             "site-dependent SU(3) links, periodic space and antiperiodic time. "
             "The links are seeded synthetic near-identity matrices, not a "
             "QCD gauge ensemble. The benchmark uses gauge strength 0.25 and "
             "positive bare masses 0.1 and 0.5.", "",
             "All three paths use **complex BiCGStab directly on D**, with "
             "independent recurrence coefficients and stopping masks for every "
             "RHS. Mixed paths use FP32 corrections (inner tolerance 0.002), "
             "FP64 accumulation, and full CPU FP64 true-residual checks. The "
             "CPU FP64 baseline uses the final tolerance directly. Every RHS "
             "must have true relative residual ≤ 1e-10. No normal equations "
             "or CG are used in these Wilson solves.", "",
             "## Complete solve runtime", "",
             "Times include all BiCGStab vector operations and reductions, "
             "operator applications, Metal encoding/submission/waits, precision "
             "conversion, FP64 solution updates, and full CPU FP64 residual "
             "checks. Reused buffer allocation, matrix upload, and compilation "
             "are excluded from warm solve time and itemized separately.", "",
             "| Mass | Lattice | RHS | CPU FP64 (ms) | CPU mixed (ms) | Metal mixed (ms) | Speedup vs faster CPU | Metal ms/RHS | RHS/s |",
             "|---:|---|---:|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        d, f, g = (med(key, b) for b in backends)
        lines.append(f"| {m:g} | {l}⁴ | {n} | {d:.2f} | {f:.2f} | {g:.2f} | {min(d,f)/g:.2f}× | {g/n:.3f} | {1000*n/g:.1f} |")
    lines += ["", "## Convergence", "",
              "Iteration counts sum the maximum per-RHS BiCGStab count at "
              "each correction step. Each iteration needs up to two Wilson "
              "applications. The GPU executes commands in groups of four "
              "iterations and includes padded work after RHS converge; the "
              "CPU can stop immediately, including after the first half of "
              "an iteration. Padded work is included in the Metal timings.", "",
              "| Mass | RHS | CPU FP64 corrections / iterations | CPU mixed corrections / iterations | Metal corrections / iterations / dispatched | Worst Metal true residual |",
              "|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        counts = []
        for b in backends:
            counts.append(" / ".join(integer_range(key, b, field) for field in
                          (("outer_steps", "inner_max_sum", "dispatched_steps") if b == "metal_mixed" else ("outer_steps", "inner_max_sum"))))
        worst = max(float(r["max_true_residual"]) for r in groups[key + ("metal_mixed",)])
        lines.append(f"| {m:g} | {n} | {counts[0]} | {counts[1]} | {counts[2]} | {worst:.3g} |")
    trace_path = Path(str(args.csv)+".trace.csv")
    with trace_path.open() as f:
        traces = list(csv.DictReader(f))
    example = min(cases, key=lambda k: (k[0], -k[1], -k[2]))
    steps = defaultdict(dict)
    for r in traces:
        if (float(r["mass"]), int(r["extent"]), int(r["rhs"])) == example and r["trial"] == "0" and r["backend"] != "cpu_fp64":
            steps[int(r["outer_step"])][r["backend"]] = float(r["max_true_residual"])
    lines += ["", f"True-residual history: mass {example[0]:g}, {example[1]}⁴, "
              f"{example[2]} RHS, first measured trial.", "",
              "| Correction | CPU mixed | Metal mixed |", "|---:|---:|---:|"]
    for step, vals in sorted(steps.items()):
        lines.append(f"| {step} | {vals.get('cpu_mixed',float('nan')):.3g} | {vals.get('metal_mixed',float('nan')):.3g} |")
    lines += ["", "## Runtime components and variation", "",
              "Component medians are independent and need not sum to the "
              "median total. Setup is the original pre-warmup allocation and "
              "upload cost. It is not a measured cold-start total.", "",
              "| Mass | RHS | CPU mixed min–max (ms) | Metal total min–max (ms) | FP32 inner (ms) | CPU FP64 checks (ms) | FP64 update (ms) | Metal setup (ms) |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        times = [float(r["total_ms"]) for r in groups[key + ("metal_mixed",)]]
        cpu_times = [float(r["total_ms"]) for r in groups[key + ("cpu_mixed",)]]
        vals = [med(key, "metal_mixed", f) for f in ("inner_ms", "residual_ms", "update_ms", "setup_ms")]
        lines.append(f"| {m:g} | {n} | {min(cpu_times):.2f}–{max(cpu_times):.2f} | {min(times):.2f}–{max(times):.2f} | {vals[0]:.2f} | {vals[1]:.2f} | {vals[2]:.2f} | {vals[3]:.2f} |")
    lib = statistics.median(float(r["library_ms"]) for r in data)
    lines += ["", f"Shared Metal library initialization was {lib:.2f} ms, "
              "also outside warm solve time. CPU setup costs and all per-trial "
              "measurements are preserved in the CSV. This is a local run with "
              "ordinary system activity, not a thermally or scheduler-isolated "
              "measurement; CPU timings in particular show variation.", "", "## RHS batching", ""]
    for m in sorted({m for m, l, n in cases}):
        subset = [k for k in cases if k[0] == m]
        best = min(subset, key=lambda k: med(k,"metal_mixed")/k[2])
        lines.append(f"For mass {m:g}, the highest measured Metal throughput is "
                     f"{1000*best[2]/med(best,'metal_mixed'):.1f} RHS/s at {best[2]} "
                     "RHS. Compare this with the timing spread above before "
                     "interpreting a small difference between batch sizes.")
        lines.append("")
    lines += ["These are independent batched solves, not block BiCGStab or "
              "recycling of a common Krylov basis. Additional RHS amortize "
              "control overhead and improve utilization, but increase vector "
              "traffic and FP64 CPU work. They do not guarantee higher "
              "throughput. The earlier SPD throughput numbers do not predict "
              "Wilson throughput or iteration counts.", "",
              "## Validation and scope", "",
              f"All {len(data)} measured complete solves met the per-RHS "
              "true-residual target. Every component was checked on the CPU "
              "in FP64 against the original FP64 Wilson operator. Known-solution "
              f"relative error was at most {max(float(r['solution_error']) for r in data):.3g}. "
              "The separate nine-system self-test checks gamma anticommutation, "
              "SU(3) unitarity and determinant, gamma5 Hermiticity, non-Hermiticity, "
              "direct Wilson versus encoded CPU/Metal products, free-field "
              "Fourier symbol, odd RHS, and a zero RHS. CPU and Metal reductions "
              "have different orders, so BiCGStab iteration counts can differ.", "",
              "No clover term, even/odd Schur complement, multigrid, domain "
              "decomposition, physical negative bare mass, or production "
              "gauge ensemble is included. No Chroma/QDP binary was used as a "
              "numerical oracle. Operator normalization and the identity/gamma "
              "Kronecker decomposition follow the inspected local source, "
              "with an explicitly documented Euclidean spin basis. This is "
              "closer to mgproton.cc's workload, not a production substitute.", "",
              "See [WILSON.md](WILSON.md) for equations, source references, "
              "implementation details, and reproduction.", "",
              f"Raw data: [{args.csv.as_posix()}]({args.csv.as_posix()}), "
              f"[convergence traces]({trace_path.as_posix()}).", ""]
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output}: {len(cases)} cases, {len(data)} solves")


if __name__ == "__main__":
    main()
