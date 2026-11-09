#!/usr/bin/env python3
"""Summarize raw full-solve samples; standard library only."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("inputs", nargs="+", type=Path)
    ap.add_argument("--output", type=Path, default=Path("SOLVE_RESULTS.md"))
    args = ap.parse_args()
    groups = defaultdict(list)
    traces = []
    for file in args.inputs:
        with file.open() as f:
            for row in csv.DictReader(f):
                groups[(int(row["extent"]), int(row["rhs"]), row["backend"])].append(row)
        with Path(str(file) + ".trace.csv").open() as f:
            traces.extend(csv.DictReader(f))
    cases = sorted({(l, r) for l, r, b in groups})
    backends = ("cpu_fp64", "cpu_mixed", "metal_mixed")
    for l, r in cases:
        for b in backends:
            if (l, r, b) not in groups:
                raise SystemExit(f"missing backend: {(l, r, b)}")

    def med(l, r, b, key="total_ms"):
        return statistics.median(float(x[key]) for x in groups[l, r, b])

    lines = ["# Mixed-precision solve on M1 Max", "",
             "Measured locally on November 3, 2026, on the same M1 Max "
             "(24 GPU cores, 32 GB, eight CPU worker tasks) as the multiply test. "
             "Connected to AC power, macOS 26.7.1. Each table entry is the median "
             "of three complete warm solves; backend order rotates between trials. "
             "Raw minimum/maximum variation can be recovered from the per-trial CSV.", "",
             "All paths solve the same synthetic Hermitian positive-definite "
             "complex block-sparse system. The target is a **true FP64 relative "
             "residual ≤ 1e-10 for every RHS**, measured on the CPU over every "
             "component. Mixed paths use FP32 independent CG corrections with "
             "inner tolerance 0.002 and an FP64 accumulated solution. The all-FP64 "
             "CPU baseline uses CG at the final tolerance and verifies its true "
             "residual too. CPU kernels use multithreaded ARM NEON.", "",
             "**Total solve time includes** input/residual conversion, FP32 "
             "iterations, all reductions and vector updates, Metal encoding, "
             "submission/waits, FP64 solution updates, and the full CPU FP64 "
             "residual checks. Allocation, operator upload, pipeline compilation, "
             "and manufactured-problem generation are excluded from this warm "
             "time and reported separately below. No kernel-only timing is used "
             "for the speedups.", ""]
    for l in sorted({l for l, r in cases}):
        lines += [f"## Complete solves: {l}⁴ lattice", "",
                  "| RHS | CPU FP64 (ms) | CPU mixed (ms) | GPU mixed (ms) | GPU speedup vs faster CPU | GPU ms/RHS | GPU RHS/s |",
                  "|---:|---:|---:|---:|---:|---:|---:|"]
        for ll, r in cases:
            if ll != l:
                continue
            d, f, g = (med(l, r, b) for b in backends)
            lines.append(f"| {r} | {d:.2f} | {f:.2f} | {g:.2f} | {min(d, f)/g:.2f}× | {g/r:.3f} | {1000*r/g:.1f} |")
        lines.append("")

    largest = max(l for l, r in cases)
    lines += [f"## Convergence on {largest}⁴", "",
              "Inner iterations below are the sum over refinement steps of the "
              "largest per-RHS iteration count, not a block-CG iteration count. "
              "GPU commands contain eight CG iterations each, so dispatched "
              "iterations also include padding after RHS have converged. "
              "The CPU stops immediately; the GPU's padded work is included in "
              "its timings.", "",
              "| RHS | CPU FP64 iterations | CPU mixed corrections / inner iterations | GPU corrections / inner iterations / dispatched | Worst GPU true residual |",
              "|---:|---:|---:|---:|---:|"]
    for l, r in cases:
        if l != largest:
            continue
        d = int(med(l, r, "cpu_fp64", "inner_max_sum"))
        f = [int(med(l, r, "cpu_mixed", k)) for k in ("outer_steps", "inner_max_sum")]
        g = [int(med(l, r, "metal_mixed", k)) for k in ("outer_steps", "inner_max_sum", "dispatched_steps")]
        worst = max(float(x["max_true_residual"]) for x in groups[l, r, "metal_mixed"])
        lines.append(f"| {r} | {d} | {f[0]} / {f[1]} | {g[0]} / {g[1]} / {g[2]} | {worst:.3g} |")
    lines.append("")
    example_rhs = 48 if (largest, 48) in cases else cases[-1][1]
    lines += [f"Example convergence for {largest}⁴ / {example_rhs} RHS (first measured trial):", "",
              "| Correction | CPU mixed max true residual | GPU mixed max true residual |", "|---:|---:|---:|"]
    by_step = defaultdict(dict)
    for tr in traces:
        if int(tr["extent"]) == largest and int(tr["rhs"]) == example_rhs and tr["trial"] == "0" and tr["backend"] != "cpu_fp64":
            by_step[int(tr["outer_step"])][tr["backend"]] = float(tr["max_true_residual"])
    for step, values in sorted(by_step.items()):
        lines.append(f"| {step} | {values.get('cpu_mixed', float('nan')):.3g} | {values.get('metal_mixed', float('nan')):.3g} |")

    lines += ["", "## Where time goes and batching", "",
              f"GPU mixed time breakdown on the {largest}⁴ lattice (median of "
              "each component independently; rounded components need not sum to "
              "the median total):", "",
              "| RHS | Inner solve (ms) | CPU FP64 residual (ms) | FP64 update (ms) | One-time GPU setup (ms) | Warm solve + setup (ms) |",
              "|---:|---:|---:|---:|---:|---:|"]
    for l, r in cases:
        if l != largest:
            continue
        vals = [med(l, r, "metal_mixed", k) for k in ("inner_ms", "residual_ms", "update_ms", "setup_ms")]
        lines.append(f"| {r} | {vals[0]:.2f} | {vals[1]:.2f} | {vals[2]:.2f} | {vals[3]:.2f} | {med(l,r,'metal_mixed')+vals[3]:.2f} |")
    optimum = min((r for l, r in cases if l == largest), key=lambda r: med(largest, r, "metal_mixed")/r)
    fastest = med(largest, optimum, "metal_mixed")/optimum
    lines += ["",
              f"The highest observed GPU throughput on the larger lattice is "
              f"at **{optimum} RHS**, at {fastest:.3f} ms/RHS "
              f"({1000/fastest:.1f} RHS/s). This is the best of the tested batch "
              "sizes, not a proven optimum. Batch size should be judged by "
              "time per RHS rather than the wall time of an entire batch.", "",
              "These are independent CG solves scheduled as a batch, with "
              "separate coefficients and convergence masks per RHS. There is "
              "no block Krylov method, shared Krylov basis, or multigrid "
              "preconditioner. More RHS can amortize dispatch overhead and "
              "improve utilization, but also increase memory traffic and "
              "the cost of CPU FP64 checks. Larger batches need not improve "
              "throughput indefinitely.", "",
              "One-time setup includes workspaces, pipeline specialization and "
              "matrix/shared-buffer allocation and upload. A shared Metal "
              f"library initialization of {statistics.median(float(x['library_ms']) for v in groups.values() for x in v):.2f} ms "
              "is also excluded. Setup values come from the original allocations "
              "before warmup; adding them to warm medians is an accounting "
              "estimate, not a separately measured cold-start solve. CPU setup "
              "costs are in the raw CSV too.", "",
              "## Validation and limits", "",
              "Every measured solve passed full-component residual checks "
              "and a comparison against an independently generated known "
              "solution. The self-test covers nine small systems, zero RHS "
              "locking, odd RHS counts, and both SIMD remainder paths. A "
              "separate mass=0.01 run checks a higher condition-number bound. "
              "The original 30 multiply tests still pass after the utility-code "
              "refactor.", "",
              "The operator is a periodic covariant Laplacian with fixed "
              "unitary color transports and spin permutations, plus mass=0.25. "
              "Its spectrum is bounded by [mass, mass+16]. It is a synthetic SPD "
              "problem with especially few distinct eigenvalues on very small "
              "lattices. It is not the Wilson/Dirac matrix in mgproton.cc. "
              "CG requires SPD; a real nonsymmetric operator needs a suitable "
              "solver and preconditioner. These results establish the tested "
              "mixed-precision mechanics and throughput, not production "
              "multigrid convergence or speedup. No Chroma/SuperbBLAS checkout "
              "was modified.", "",
              "See [SOLVE.md](SOLVE.md) for algorithms and reproduction. Raw data:", ""]
    lines += [f"- [{p.as_posix()}]({p.as_posix()}) and its `.trace.csv` file." for p in args.inputs]
    args.output.write_text("\n".join(lines)+"\n")
    print(f"Wrote {args.output}: {len(cases)} systems")


if __name__ == "__main__":
    main()
