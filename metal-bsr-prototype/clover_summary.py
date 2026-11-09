#!/usr/bin/env python3
"""Summarize c_sw=1 odd-Schur solves and the elimination precision change."""
import argparse
import csv
import statistics
from collections import defaultdict
from pathlib import Path


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv", type=Path)
    ap.add_argument("--output", type=Path, default=Path("CLOVER_RESULTS.md"))
    args = ap.parse_args()
    with args.csv.open() as f:
        data = list(csv.DictReader(f))
    if not data:
        raise SystemExit("No measurements")
    groups = defaultdict(list)
    for row in data:
        if float(row["max_true_residual"]) > float(row["tolerance"]):
            raise SystemExit("Full clover-Wilson residual target missed")
        key = (float(row["mass"]), int(row["extent"]), int(row["rhs"]))
        groups[key + (row["backend"],)].append(row)
    cases = sorted({k[:3] for k in groups})
    backends = ("cpu_fp64", "cpu_mixed_host64", "cpu_mixed_fp32", "metal_mixed_host64", "metal_mixed_fp32")
    for key in cases:
        for b in backends:
            if len(groups[key + (b,)]) != 3:
                raise SystemExit(f"Expected three trials for {key} / {b}")

    def med(key, b, field="total_ms"):
        return statistics.median(float(r[field]) for r in groups[key + (b,)])

    def count(key, b):
        values = []
        for field in ("outer_steps", "inner_max_sum", "dispatched_steps"):
            v = sorted({int(r[field]) for r in groups[key + (b,)]})
            values.append(str(v[0]) if len(v) == 1 else f"{v[0]}–{v[-1]}")
        return " / ".join(values)

    lines = ["# Clover odd-Schur / Metal FP32 elimination results", "",
             "Measured locally on November 4, 2026, on the M1 Max (24 GPU cores, "
             "32 GB, eight CPU worker tasks), macOS 26.7.1, connected to AC power. "
             "Three complete solves from zero are measured per method after "
             "one full warmup, with rotating method order. All methods in a "
             "case share identical gauge links, clover blocks, and RHS.", "",
             "The operator uses **c_sw=1.0**, QDP's gamma matrices with diagonal "
             "gamma5, and **two Hermitian 6×6 chiral blocks at every site**. "
             "The retained odd system is exactly `S_o = A_oo - H_oe A_ee^-1 H_eo`, "
             "without scalar mass division or Schur normalization. H includes "
             "the Wilson `-1/2 Dslash` factor. The links are seeded synthetic "
             "SU(3) matrices at strength 0.25, with periodic space and "
             "antiperiodic fermionic time boundaries.", "",
             "Chroma's normalization is explicit: `F_AH=(Q-Q†)/8` from four "
             "oriented leaves and `C=-(c_sw/2) sum_(mu<nu) gamma_mu gamma_nu F_AH`. "
             "This agrees with the inspected Chroma packed-component formulas. "
             "The convention, local inverses, and full clover solves were "
             "validated before implementing the Metal elimination change.", "",
             "Every solve is accepted only when **every RHS has full CPU FP64 "
             "true relative residual ≤ 1e-10 after even-site reconstruction**. "
             "FP64 accumulated solutions and checks use the original FP64 "
             "clover-Wilson coefficients. Mixed inner BiCGStab tolerance is 0.002.", "",
             "## Complete solve runtime", "",
             "All times include input conversion, elimination, Schur BiCGStab, "
             "reductions and vector updates, GPU submission/waits, reconstruction, "
             "FP64 accumulation, and full FP64 residual checks. Allocation, "
             "clover construction/inversion, and compilation are outside warm "
             "time and recorded separately. CPU best is the faster of the three "
             "measured CPU methods, not a tuned production Dslash baseline.", "",
             "| Mass | Lattice | RHS | CPU best (ms) | CPU mixed FP32 elim. (ms) | Metal CPU-FP64 elim. (ms) | Metal device-FP32 elim. (ms) | Move gain | CPU best / Metal new | New ms/RHS |",
             "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        best = min(med(key, b) for b in backends[:3])
        cpu, old, new = (med(key, b) for b in ("cpu_mixed_fp32", "metal_mixed_host64", "metal_mixed_fp32"))
        lines.append(f"| {m:g} | {l}⁴ | {n} | {best:.2f} | {cpu:.2f} | {old:.2f} | {new:.2f} | {old/new:.2f}× | {best/new:.2f}× | {new/n:.3f} |")
    lines += ["", "## All methods and convergence", "",
              "Names ending in `host64` use CPU FP64 RHS reduction and "
              "reconstruction around FP32 inner solves. `cpu_mixed_fp32` "
              "performs all elimination math in CPU FP32. `metal_mixed_fp32` "
              "does that math, packing/scattering, and inner initialization "
              "on Metal. All paths retain CPU FP64 true-residual checks. "
              "Counts are `corrections / summed maximum inner iterations / "
              "dispatched iterations`; Metal pads to four-iteration chunks.", "",
              "| Mass | RHS | Method | Total (ms) | Corrections / inner / dispatched | Worst true residual | Known-solution error |",
              "|---:|---:|---|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        for b in backends:
            rows = groups[key + (b,)]
            residual = max(float(r["max_true_residual"]) for r in rows)
            error = max(float(r["solution_error"]) for r in rows)
            lines.append(f"| {m:g} | {n} | {b} | {med(key,b):.2f} | {count(key,b)} | {residual:.3g} | {error:.3g} |")
    lines += ["", "## Time breakdown of the Metal change", "",
              "Component medians are independent and may not sum to total. "
              "The inner column subtracts preparation/reconstruction subfields "
              "from raw `inner_ms`, avoiding double counting.", "",
              "| Mass | RHS | Elimination | Schur inner (ms) | Preparation (ms) | Reconstruction (ms) | CPU FP64 residual (ms) | FP64 update (ms) |",
              "|---:|---:|---|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        for b, label in (("metal_mixed_host64", "CPU FP64"), ("metal_mixed_fp32", "Metal FP32")):
            rows = groups[key + (b,)]
            pure = statistics.median(float(r["inner_ms"])-float(r["prepare_ms"])-float(r["reconstruct_ms"]) for r in rows)
            vals = [med(key,b,f) for f in ("prepare_ms", "reconstruct_ms", "residual_ms", "update_ms")]
            lines.append(f"| {m:g} | {n} | {label} | {pure:.2f} | " + " | ".join(f"{v:.2f}" for v in vals) + " |")
    lines += ["", "Metal preparation still includes a CPU FP64-to-FP32 cast of "
              "the full residual into shared memory, because Metal has no native "
              "FP64 arithmetic. The reduced odd RHS stays on the GPU; it is "
              "not downloaded to initialize BiCGStab. The FP32 full correction "
              "is read only after GPU reconstruction completes, then accumulated "
              "in CPU FP64. End-to-end timing includes these boundaries.", ""]
    trace_path = Path(str(args.csv)+".trace.csv")
    with trace_path.open() as f:
        traces = list(csv.DictReader(f))
    example = min(cases,key=lambda k:(k[0],-k[1],-k[2]))
    history = defaultdict(dict)
    for tr in traces:
        if (float(tr["mass"]),int(tr["extent"]),int(tr["rhs"])) == example and tr["trial"] == "0" and tr["backend"].startswith("metal"):
            history[int(tr["outer_step"])][tr["backend"]] = float(tr["max_true_residual"])
    lines += [f"## Full FP64 residual history: mass {example[0]:g}, {example[2]} RHS", "",
              "| Correction | Metal + CPU FP64 elimination | Metal FP32 elimination |", "|---:|---:|---:|"]
    for step, vals in sorted(history.items()):
        lines.append(f"| {step} | {vals.get('metal_mixed_host64',float('nan')):.3g} | {vals.get('metal_mixed_fp32',float('nan')):.3g} |")
    lines += ["", "## Batching and timing spread", "",
              "| Mass | RHS | New Metal min–max (ms) | New Metal RHS/s | Workspace setup (ms) | Problem setup (ms) |",
              "|---:|---:|---:|---:|---:|---:|"]
    for key in cases:
        m, l, n = key
        times = [float(r["total_ms"]) for r in groups[key + ("metal_mixed_fp32",)]]
        lines.append(f"| {m:g} | {n} | {min(times):.2f}–{max(times):.2f} | {1000*n/statistics.median(times):.1f} | "
                     f"{med(key,'metal_mixed_fp32','workspace_setup_ms'):.2f} | {med(key,'metal_mixed_fp32','problem_setup_ms'):.2f} |")
    lib = statistics.median(float(r["library_ms"]) for r in data)
    lines += ["", f"Metal library initialization was {lib:.2f} ms. Problem setup "
              "includes links, clover construction, FP64 inverses, and manufactured "
              "RHS generation; workspace setup includes per-method buffers and "
              "pipelines. Shared resources are charged to each method separately, "
              "so setup columns should not be summed across methods.", "",
              "Compare time per RHS and the timing spread when choosing a batch. "
              "These are independent batched BiCGStab solves, not a shared or "
              "block Krylov solver. The run used ordinary system activity and "
              "was not thermally or scheduler isolated.", "", "## Validation and limits", "",
              f"All {len(data)} measured solves passed full FP64 residual and "
              "known-solution checks. Stage 1 validated Chroma normalization "
              "and packed formulas, Hermitian chiral blocks, independent local "
              "inverses, gauge covariance, direct full clover and odd Schur "
              "products, and complete solves before device elimination changed. "
              "The final self-test adds direct device RHS/reconstruction checks "
              "and all five solver paths on eight small cases, including odd "
              "RHS, a zero RHS, free field, and c_sw=0. An analytic quantized "
              "Abelian field independently fixes the plaquette normalization "
              "and orientation.", "",
              "The FP64 inverses are precomputed once and rounded for FP32 "
              "application; there is no GPU factorization claim. Clover/local "
              "matrices are reused across RHS. Field strengths are built once "
              "from gauge links; adding clover does not automatically reduce "
              "Dslash link traffic. The data remain synthetic, isotropic, and "
              "positive bare mass with no multigrid, production ensemble, or "
              "negative-critical-mass validation. Generic kernels are used on "
              "both CPU and GPU. The historical scalar-diagonal Wilson results "
              "are not an apples-to-apples baseline for the new clover operator.", "",
              "See [CLOVER.md](CLOVER.md) for exact equations and reproduction.", "",
              f"Raw data: [{args.csv.as_posix()}]({args.csv.as_posix()}), "
              f"[convergence traces]({trace_path.as_posix()}).", ""]
    args.output.write_text("\n".join(lines))
    print(f"Wrote {args.output}: {len(cases)} cases / {len(data)} solves")


if __name__ == "__main__":
    main()
