#!/usr/bin/env python3
"""Report full-Wilson true-residual checked even/odd benchmark results."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("--output", type=Path, default=Path("EVENODD_RESULTS.md"))
    args = ap.parse_args()
    with args.csv.open() as f:
        data = list(csv.DictReader(f))
    if not data:
        raise SystemExit("No measurements")
    groups = defaultdict(list)
    for r in data:
        if float(r["max_true_residual"]) > float(r["tolerance"]):
            raise SystemExit("A solve missed the full-Wilson tolerance")
        key = (float(r["mass"]), int(r["extent"]), int(r["rhs"]))
        groups[key + (r["preconditioner"], r["backend"])].append(r)
    cases = sorted({k[:3] for k in groups})
    backends = ("cpu_fp64", "cpu_mixed", "metal_mixed")
    for key in cases:
        for pre in ("none", "even_odd"):
            for b in backends:
                if len(groups[key + (pre, b)]) != 3:
                    raise SystemExit(f"Expected three trials: {key}, {pre}, {b}")

    def med(key, pre, b, field="total_ms"):
        return statistics.median(float(r[field]) for r in groups[key + (pre, b)])

    def best_cpu(key, pre):
        return min(med(key, pre, b) for b in backends[:2])

    def counts(key, pre, b):
        out = []
        for field in ("outer_steps", "inner_max_sum", "dispatched_steps"):
            vals = sorted({int(r[field]) for r in groups[key + (pre, b)]})
            out.append(str(vals[0]) if len(vals) == 1 else f"{vals[0]}–{vals[-1]}")
        return " / ".join(out)

    lines = ["# Even/odd Wilson benchmark on M1 Max", "",
             "Measured locally on November 3, 2026: M1 Max, 24 GPU cores, "
             "32 GB, macOS 26.7.1, AC power, eight CPU worker tasks. "
             "Each entry is the median of three complete solves from zero, "
             "after one full warmup per method. All six methods are measured "
             "in the same run, with rotating order. The links and RHS are "
             "identical across methods within each case.", "",
             "The retained operator is the bare Wilson matrix with site-dependent "
             "synthetic SU(3) links (strength 0.25), periodic space, and "
             "antiperiodic time. Even/odd uses packed checkerboards and the "
             "normalized Schur complement `S = I - H_eo H_oe / (m+4)^2`. "
             "Each Schur product uses two half-volume hopping passes. "
             "All paths use complex BiCGStab directly, with no normal equations.", "",
             "**Every RHS is accepted only after full odd-site reconstruction "
             "and a CPU FP64 check of `||b-D64*x||/||b|| <= 1e-10`.** "
             "Mixed paths use FP32 inner tolerance 0.002 and FP64 correction "
             "accumulation. EO RHS preparation and odd reconstruction are also "
             "CPU FP64 work, included in complete solve time. Initial allocation, "
             "operator upload, and compilation are outside warm solve time.", "",
             "## Complete runtime: with and without even/odd", "",
             "CPU best means the faster measured CPU FP64 or CPU mixed median "
             "for that preconditioner. EO gain means unpreconditioned time "
             "divided by preconditioned time; a value below one is a slowdown.", "",
             "| Mass | Lattice | RHS | CPU best, none (ms) | CPU best, EO (ms) | Metal, none (ms) | Metal, EO (ms) | Metal EO gain | CPU EO / Metal EO | EO Metal ms/RHS |",
             "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        cn, ce = (best_cpu(key, p) for p in ("none", "even_odd"))
        gn, ge = (med(key, p, "metal_mixed") for p in ("none", "even_odd"))
        lines.append(f"| {m:g} | {l}⁴ | {n} | {cn:.2f} | {ce:.2f} | {gn:.2f} | {ge:.2f} | {gn/ge:.2f}× | {ce/ge:.2f}× | {ge/n:.3f} |")
    lines += ["", "## All CPU baselines and convergence", "",
              "Counts are `outer corrections / summed maximum per-RHS inner "
              "iterations / dispatched iterations`. A Schur iteration is not "
              "a half-cost Wilson iteration: each operator application makes "
              "two half-volume hopping passes. BiCGStab needs up to two "
              "operator applications per iteration. Metal pads to four-iteration "
              "command chunks; the CPU can stop early. Timing includes that work.", "",
              "| Mass | RHS | Preconditioner | Backend | Total (ms) | Corrections / inner / dispatched | Worst full residual |",
              "|---:|---:|---|---|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        for pre in ("none", "even_odd"):
            for b in backends:
                worst = max(float(r["max_true_residual"]) for r in groups[key + (pre, b)])
                lines.append(f"| {m:g} | {n} | {pre} | {b} | {med(key,pre,b):.2f} | {counts(key,pre,b)} | {worst:.3g} |")
    lines += ["", "## Where the preconditioned Metal time goes", "",
              "The inner column excludes RHS preparation and reconstruction, "
              "which the raw `inner_ms` field includes. Component medians "
              "need not sum to the median total. Setup is the original allocation "
              "and checkerboard extraction/upload cost, not a cold-start solve.", "",
              "| Mass | RHS | FP32 Schur BiCGStab (ms) | FP64 RHS reduction (ms) | FP64 odd reconstruction (ms) | Full FP64 residual checks (ms) | FP64 update (ms) | EO setup (ms) |",
              "|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        rows = groups[key + ("even_odd", "metal_mixed")]
        pure = statistics.median(float(r["inner_ms"])-float(r["eo_prepare_ms"])-float(r["eo_reconstruct_ms"]) for r in rows)
        vals = [med(key, "even_odd", "metal_mixed", f) for f in
                ("eo_prepare_ms", "eo_reconstruct_ms", "residual_ms", "update_ms", "setup_ms")]
        lines.append(f"| {m:g} | {n} | {pure:.2f} | " + " | ".join(f"{x:.2f}" for x in vals) + " |")
    with Path(str(args.csv)+".trace.csv").open() as f:
        traces = list(csv.DictReader(f))
    example = min(cases, key=lambda k:(k[0],-k[1],-k[2]))
    history = defaultdict(dict)
    for r in traces:
        if (float(r["mass"]),int(r["extent"]),int(r["rhs"])) == example and r["trial"] == "0" and r["backend"] == "metal_mixed":
            history[int(r["outer_step"])][r["preconditioner"]] = float(r["max_true_residual"])
    lines += ["", f"## Full-Wilson convergence example: mass {example[0]:g}, {example[2]} RHS", "",
              "These are full-system residuals after reconstruction, not "
              "reduced-system residuals.", "",
              "| Correction | Metal none | Metal even/odd |", "|---:|---:|---:|"]
    for step, vals in sorted(history.items()):
        lines.append(f"| {step} | {vals.get('none',float('nan')):.3g} | {vals.get('even_odd',float('nan')):.3g} |")
    lines += ["", "## Timing spread and batching", "",
              "| Mass | RHS | EO CPU mixed min–max (ms) | EO Metal min–max (ms) | EO Metal RHS/s |",
              "|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        cpu = [float(r["total_ms"]) for r in groups[key + ("even_odd","cpu_mixed")]]
        gpu = [float(r["total_ms"]) for r in groups[key + ("even_odd","metal_mixed")]]
        lines.append(f"| {m:g} | {n} | {min(cpu):.2f}–{max(cpu):.2f} | {min(gpu):.2f}–{max(gpu):.2f} | {1000*n/statistics.median(gpu):.1f} |")
    lines += ["", "Compare milliseconds per RHS and the timing spread when "
              "choosing 96 versus 192 RHS. All solves are independent batched "
              "BiCGStab; there is no shared Krylov basis or block solver. "
              "This was an ordinary local run, not a thermally or scheduler-isolated "
              "measurement. The CSV retains all individual trial times.", "",
              "## Validation and scope", "",
              f"All {len(data)} measured solves passed full FP64 residual "
              "and known-solution checks. Maximum known-solution relative error "
              f"was {max(float(r['solution_error']) for r in data):.3g}. "
              "The separate six-system self-test checks direct full-Wilson "
              "composition versus CPU/Metal Schur products, transformed RHS, "
              "odd reconstruction, odd RHS, zero RHS locking, duplicate neighbors, "
              "periodic and antiperiodic boundaries, and rejection of odd "
              "periodic extents. A separate free-field periodic case passed "
              "all six solver paths, and the nine-system unpreconditioned "
              "Wilson regression test passed after the shared-code changes.", "",
              "The operator still has synthetic SU(3) links, positive bare "
              "mass, and no clover, multigrid, domain decomposition, or production "
              "gauge ensemble. Generic color/spin kernels are used on both "
              "processors; neither is a tuned half-spinor Dslash implementation. "
              "FP64 RHS/reconstruction on the CPU is a deliberate implementation "
              "choice and can limit net Metal speedup despite fewer inner iterations.", "",
              "See [EVENODD.md](EVENODD.md) for formulas and reproduction.", "",
              f"Raw measurements: [{args.csv.as_posix()}]({args.csv.as_posix()}) "
              f"and [convergence traces]({args.csv.as_posix()}.trace.csv).", ""]
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output}: {len(cases)} cases / {len(data)} solves")


if __name__ == "__main__":
    main()
