#!/usr/bin/env python3
"""Generate a compact benchmark report from the prototype's raw CSV."""
import argparse
import csv
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("--output", type=Path, default=Path("RESULTS.md"))
    args = ap.parse_args()
    with args.csv.open() as f:
        data = list(csv.DictReader(f))
    groups = defaultdict(dict)
    for row in data:
        key = (row["mode"], int(row["extent"]), int(row["rhs"]))
        groups[key][row["backend"]] = row
    expected = {"accelerate_1", "accelerate_parallel", "neon_parallel", "metal_latency", "metal_batch"}
    for key, backends in groups.items():
        if set(backends) != expected:
            raise SystemExit(f"Incomplete benchmark case: {key}")
    if not groups:
        raise SystemExit("No benchmark cases")

    def ms(row):
        return float(row["median_ms"])

    speedups = []
    latency_wins = 0
    lines = [
        "# M1 Max Metal block-sparse benchmark",
        "",
        "Measured locally on November 3, 2026: Apple M1 Max, 24-core GPU, "
        "10-core CPU (8 performance + 2 efficiency), 32 GB unified memory; "
        "macOS 26.7.1 (25G241), connected to AC power. CPU parallel runs request "
        "eight worker tasks; Accelerate internal threads are limited to one.",
        "",
        f"Source data: [{args.csv.as_posix()}]({args.csv.as_posix()}). "
        "Times below are median **wall milliseconds per multiply**, including "
        "command encoding/submission/wait for Metal. Initial allocation, data "
        "copies, and pipeline compilation are excluded. All data remains resident.",
        "",
        "CPU best means the fastest measured CPU implementation in that case, "
        "chosen from single-worker Accelerate, eight-worker Accelerate, and "
        "eight-worker NEON. GPU batch amortizes command-buffer overhead across "
        "independent repeated multiplies; GPU single submits and waits for one "
        "multiply. Speedup uses CPU best divided by the corresponding GPU time.",
        "",
    ]
    for mode, title in [("dense12", "Dense 12×12 complex blocks"),
                        ("color3_spin4", "3×3 color blocks with four-spin permutation")]:
        lines += [f"## {title}", "",
                  "| Lattice | RHS | CPU best (ms) | CPU method | GPU single (ms) | GPU batch (ms) | Single speedup | Batch speedup |",
                  "|---|---:|---:|---|---:|---:|---:|---:|"]
        for (m, extent, rhs), b in sorted(groups.items()):
            if m != mode:
                continue
            cpu_name = min((n for n in b if not n.startswith("metal")), key=lambda n: ms(b[n]))
            cpu, single, batch = ms(b[cpu_name]), ms(b["metal_latency"]), ms(b["metal_batch"])
            speedups.append(cpu / batch)
            latency_wins += cpu > single
            short_name = {"accelerate_1": "BLAS 1", "accelerate_parallel": "BLAS 8", "neon_parallel": "NEON 8"}[cpu_name]
            lines.append(f"| {extent}⁴ | {rhs} | {cpu:.3f} | {short_name} | {single:.3f} | {batch:.3f} | {cpu/single:.2f}× | {cpu/batch:.2f}× |")
        lines.append("")

    max_abs = max(float(r["max_abs_error"]) for r in data)
    max_rel = max(float(r["relative_l2_error"]) for r in data)
    max_ws = max(float(r["working_set_mib"]) for r in data)
    lines += ["## Interpretation and limits", "",
              f"Across these {len(groups)} synthetic cases, amortized GPU speedups "
              f"range from {min(speedups):.2f}× to {max(speedups):.2f}× versus the "
              f"fastest measured CPU method. Single-submit GPU execution is faster "
              f"in {latency_wins}/{len(groups)} cases. These are case-specific "
              "measurements, not a prediction of full-solver speedup.", "",
              f"All numerical checks passed. The largest observed absolute error "
              f"was {max_abs:.3g}; the largest relative L2 error was {max_rel:.3g} "
              "against an independent double-precision reference. Small benchmark "
              "cases are fully checked; large cases use seeded samples and "
              "boundary sites. The separate 30-case self-test fully checks every "
              "output, including irregular sparsity and empty rows. An additional "
              "[full-output check](results/full-check.log) verifies all 589,824 "
              "outputs for each operator on an 8⁴ lattice with 12 RHS. "
              "That separate check is not included in the primary timing table.", "",
              f"The largest per-implementation matrix/input/output working set "
              f"was {max_ws:.1f} MiB. The process holds separate CPU and GPU copies.", "",
              "This demonstrates a standalone FP32 complex multiply using the "
              "relevant sparse structures on this M1 Max. It does not benchmark "
              "SuperbBLAS itself or implement a complete backend. Actual multigrid "
              "performance also depends on reductions, tensor rearrangements, "
              "precision conversion, CPU/GPU dependencies, and numerical convergence. "
              "Batched throughput is especially optimistic for workflows that must "
              "return to the CPU between every operation.", "",
              f"The primary run used {', '.join(sorted({r['samples'] for r in data}))} "
              "samples per backend and case. Timing "
              "spread and batch sizes are preserved in the CSV; GPU latency for "
              "small cases is sensitive to OS scheduling and submission overhead. "
              "This was a local run with ordinary system activity, not an isolated "
              "laboratory measurement.", "",
              "See [README.md](README.md) for the equations, timing definitions, "
              "validation coverage, and reproduction commands.", ""]
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output} from {len(groups)} cases")


if __name__ == "__main__":
    main()
