#!/usr/bin/env python3
"""Regenerate the tuning report from measured data, checking every result."""
import csv
import math
import statistics as st
from collections import defaultdict
from pathlib import Path

root = Path("results/mg-tuned")
with (root / "comparison.csv").open() as f:
    rows = list(csv.DictReader(f))
groups = defaultdict(list)
for r in rows:
    for field in ("total_ms", "max_true_residual", "solution_error"):
        assert math.isfinite(float(r[field])), (field, r)
    assert float(r["max_true_residual"]) <= float(r["tolerance"]) == 1e-10
    assert float(r["solution_error"]) <= 1e-7
    assert float(r["csw"]) == 1 and int(r["extent"]) == 8
    assert int(r["pre_smooth"]) == int(r["post_smooth"]) == 1
    assert int(r["vectors"]) == (6 if r["configuration"] == "original" else 12)
    key = (float(r["mass"]), int(r["rhs"]))
    groups[key + (r["configuration"], r["backend"])].append(r)
cases = sorted({k[:2] for k in groups})
methods = {"original": ("cpu_bicgstab", "metal_bicgstab", "cpu_fgmres_mg", "metal_fgmres_mg"),
           "tuned": ("cpu_fgmres_smooth", "metal_fgmres_smooth", "cpu_fgmres_mg", "metal_fgmres_mg")}
assert cases == [(.1, 96), (.1, 192), (.5, 96), (.5, 192)]
for key in cases:
    for config, names in methods.items():
        for name in names:
            assert len(groups[key + (config, name)]) == 3


def vals(k, config, backend, field):
    return [float(r[field]) for r in groups[k + (config, backend)]]


def med(k, config, backend, field="total_ms"):
    return st.median(vals(k, config, backend, field))


def count(k, config, backend, field):
    v = sorted({int(r[field]) for r in groups[k + (config, backend)]})
    return str(v[0]) if len(v) == 1 else f"{v[0]}–{v[-1]}"


def text_count(k, config, backend):
    return "/".join(count(k, config, backend, f) for f in
                    ("outer_steps", "inner_max_sum", "fine_schur_calls", "coarse_calls"))


lines = ["# Clover multigrid tuning on M1 Max", "",
         "Measured locally on November 4, 2026, on an Apple M1 Max, 24 GPU cores, "
         "32 GB, macOS 26.7.1. This extends only the standalone prototype. "
         "The fine operator remains the Chroma-convention odd clover Schur "
         "complement with c_sw=1, positive bare masses 0.1/0.5, synthetic "
         "gauge strength 0.25, an 8⁴ lattice, and antiperiodic time. "
         "All CPU/Metal inner arithmetic and elimination are FP32. "
         "Every full residual component, including reconstructed even sites "
         "and clover terms, is evaluated on the CPU in FP64; every RHS must "
         "meet a relative residual of 1e-10.", "",
         "## Selected settings", "",
         "The selected hierarchy doubles the local basis from 6 to **12 vectors "
         "per chirality**, retains **4×4×4×4 aggregates**, and changes the cycle's "
         "clover-block Jacobi weight from **0.65 to 1.0**. One pre- and one "
         "post-smoothing step remain. The coarse matrix grows from 192² to 384². "
         "Eight FP64 test-vector filtering sweeps at the independent setup "
         "weight 0.65 suffice: longer filtering did not reduce the measured "
         "large-batch iteration counts. FGMRES restart stays 8, the inner "
         "target stays 0.002, and full FP64 acceptance stays 1e-10.", "",
         "## Fine work and convergence", "",
         "Counts are summed maximum per-RHS Krylov iterations across outer "
         "corrections, and actual batched fine Schur applications. CPU and Metal "
         "use identical counts for these multigrid cases. One Schur application "
         "includes both half-volume hopping operations and the local clover "
         "blocks. Counts include smoothing, Arnoldi work, and residual refreshes.", "",
         "| Mass | RHS | Old → tuned Krylov iterations | Old → tuned fine S calls | Fine work reduction | Outer FP64 corrections |",
         "|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    old = med(k, "original", "metal_fgmres_mg", "fine_schur_calls")
    new = med(k, "tuned", "metal_fgmres_mg", "fine_schur_calls")
    lines.append(f"| {k[0]:g} | {k[1]} | {count(k,'original','metal_fgmres_mg','inner_max_sum')} → {count(k,'tuned','metal_fgmres_mg','inner_max_sum')} | {old:g} → {new:g} | {100*(1-new/old):.1f}% | {count(k,'tuned','metal_fgmres_mg','outer_steps')} |")

lines += ["", "## Fresh complete solve timings", "",
          "Medians of three measured solves from zero, following a full warmup "
          "for every method. Method order rotates within trials. Original and "
          "tuned configuration order alternates across mass/RHS cases; they "
          "run sequentially, never as competing benchmarks. CPU uses eight "
          "GCD workers and single-threaded Accelerate for coarse GEMM. Times "
          "include FP32 preparation/reconstruction, vector/reduction work, CPU "
          "small QR, GPU submissions/waits, FP64 accumulation, and full FP64 "
          "residual checks. Hierarchy/workspace/problem setup is excluded and "
          "reported separately. Ordinary local system activity was present.", "",
          "| Mass | RHS | CPU old MG ms | CPU tuned MG ms | Metal old MG ms | Metal tuned MG ms | CPU tuning gain | Metal tuning gain | Tuned CPU / Metal |",
          "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    co, cn, go, gn = (med(k,c,b) for c,b in (("original","cpu_fgmres_mg"),("tuned","cpu_fgmres_mg"),("original","metal_fgmres_mg"),("tuned","metal_fgmres_mg")))
    lines.append(f"| {k[0]:g} | {k[1]} | {co:.2f} | {cn:.2f} | {go:.2f} | {gn:.2f} | {co/cn:.2f}× | {go/gn:.2f}× | {cn/gn:.2f}× |")
lines += ["", "### BiCGStab and the smoother-only control", "",
          "BiCGStab is freshly measured with the original FP32 CPU/device "
          "elimination paths. It does not use the hierarchy. The smoother-only "
          "FGMRES control uses exactly the tuned two-step smoother at weight "
          "1.0, with no coarse correction. This distinguishes stronger "
          "smoothing from the coarse-grid contribution.", "",
          "| Mass | RHS | CPU BiCG ms | CPU smooth FG ms | CPU tuned MG ms | Metal BiCG ms | Metal smooth FG ms | Metal tuned MG ms |",
          "|---:|---:|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    times = [med(k,c,b) for c,b in (("original","cpu_bicgstab"),("tuned","cpu_fgmres_smooth"),("tuned","cpu_fgmres_mg"),("original","metal_bicgstab"),("tuned","metal_fgmres_smooth"),("tuned","metal_fgmres_mg"))]
    lines.append(f"| {k[0]:g} | {k[1]} | " + " | ".join(f"{t:.2f}" for t in times) + " |")
lines += ["", "A gain below 1 means the tuned two-level method is slower than "
          "that baseline. More available RHS can amortize setup only when a "
          "warm-time gain actually exists.", "",
          "| Mass | RHS | CPU BiCG / tuned MG | Metal BiCG / tuned MG | CPU smooth FG / tuned MG | Metal smooth FG / tuned MG |",
          "|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    cn, gn = (med(k,"tuned",b) for b in ("cpu_fgmres_mg","metal_fgmres_mg"))
    gains = [med(k,"original","cpu_bicgstab")/cn,med(k,"original","metal_bicgstab")/gn,
             med(k,"tuned","cpu_fgmres_smooth")/cn,med(k,"tuned","metal_fgmres_smooth")/gn]
    lines.append(f"| {k[0]:g} | {k[1]} | " + " | ".join(f"{t:.2f}×" for t in gains) + " |")

lines += ["", "## Full counts and numerical checks", "",
          "Counts below are `outer corrections / Krylov iterations / fine S "
          "applications / coarse corrections`. BiCGStab and FGMRES iteration "
          "counts represent different work; GPU BiCGStab includes dispatch "
          "padding. FGMRES refreshes its actual FP32 Schur residual on restart, "
          "and the outer solver always checks the full FP64 equation.", "",
          "| Mass | RHS | Setting | Method | Counts | Worst full FP64 residual |",
          "|---:|---:|---|---|---:|---:|"]
for k in cases:
    for config, names in methods.items():
        for b in names:
            worst = max(vals(k, config, b, "max_true_residual"))
            lines.append(f"| {k[0]:g} | {k[1]} | {config} | {b} | {text_count(k,config,b)} | {worst:.3g} |")

with (root / "comparison.trace.csv").open() as f:
    traces = list(csv.DictReader(f))
lines += ["", "Example full-system FP64 residual history, mass 0.1, 192 RHS, "
          "Metal, trial 0:", "", "| Correction | Original MG | Tuned MG | Tuned smoother FG | BiCGStab |", "|---:|---:|---:|---:|---:|"]
history = defaultdict(dict)
for r in traces:
    if float(r["mass"]) == .1 and r["rhs"] == "192" and r["trial"] == "0":
        history[int(r["outer_step"])][r["configuration"],r["backend"]] = float(r["max_true_residual"])
for i, r in sorted(history.items()):
    v = [r[c,b] for c,b in (("original","metal_fgmres_mg"),("tuned","metal_fgmres_mg"),("tuned","metal_fgmres_smooth"),("original","metal_bicgstab"))]
    lines.append(f"| {i} | " + " | ".join(f"{x:.3g}" for x in v) + " |")

lines += ["", "## Setup, throughput, and amortization", "",
          "Hierarchy setup is CPU FP64 and independent of the benchmark RHS "
          "and known solution. Each final case rebuilds it to alternate the "
          "configuration order, then reuses it for all warmups/trials/methods. "
          "The large-batch finalist runs additionally reuse one mass-specific "
          "hierarchy across 96 and 192 RHS. It can be retained for subsequent "
          "batches on the same gauge, mass, clover coefficient, and geometry.", "",
          "| Mass | RHS | Old hierarchy ms | Tuned hierarchy ms | Tuned Metal ms/RHS | Tuned Metal RHS/s | Tuned Metal FP64 residual ms | Hierarchy-only break-even batches vs Metal BiCG |",
          "|---:|---:|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    old, setup = (med(k,c,"metal_fgmres_mg","hierarchy_setup_ms") for c in ("original","tuned"))
    gn, gb = med(k,"tuned","metal_fgmres_mg"), med(k,"original","metal_bicgstab")
    be = f"{setup/(gb-gn):.1f}" if gb > gn else "No warm gain"
    lines.append(f"| {k[0]:g} | {k[1]} | {old:.2f} | {setup:.2f} | {gn/k[1]:.3f} | {1000*k[1]/gn:.1f} | {med(k,'tuned','metal_fgmres_mg','residual_ms'):.2f} | {be} |")
lines += ["", "Break-even is hierarchy setup divided by measured warm time "
          "saved per batch, without compilation/workspace setup; it is an "
          "estimate, particularly sensitive when the warm-time difference is "
          "small. CPU FP64 residual time remains inside all warm totals. "
          "Problem/workspace setup and per-trial timing components are in the raw CSV.", "",
          "| Mass | RHS | CPU old MG min–max ms | CPU tuned MG min–max ms | Metal old MG min–max ms | Metal tuned MG min–max ms |",
          "|---:|---:|---:|---:|---:|---:|"]
for k in cases:
    cpu, gpu = (vals(k,"tuned",b,"total_ms") for b in ("cpu_fgmres_mg","metal_fgmres_mg"))
    oldcpu, oldgpu = (vals(k,"original",b,"total_ms") for b in ("cpu_fgmres_mg","metal_fgmres_mg"))
    lines.append(f"| {k[0]:g} | {k[1]} | {min(oldcpu):.2f}–{max(oldcpu):.2f} | {min(cpu):.2f}–{max(cpu):.2f} | {min(oldgpu):.2f}–{max(oldgpu):.2f} | {min(gpu):.2f}–{max(gpu):.2f} |")
lines += ["", "The old CPU mass-0.1/96-RHS timing has substantial spread; "
          "its 2.56× median tuning gain should be read with that variability. "
          "Fine Schur counts and numerical convergence are reproducible. "
          "At mass 0.5 the tuned smoother-only method remains faster than MG. "
          "The small Metal MG advantage over BiCGStab at 192 RHS/mass 0.5 "
          "is within the observed timing spread and does not establish a "
          "reliable setup break-even there."]

screen = []
for filename in ("screen.csv", "screen2.csv", "screen3.csv"):
    with (Path("results/mg-screen") / filename).open() as f:
        screen += list(csv.DictReader(f))
assert all(float(r["max_true_residual"]) <= 1e-10 for r in screen)
configs = {r["configuration"]: r for r in screen if r["backend"] == "metal_fgmres_mg"}
names = ["original", "pre_only", "post_only", "filter24", "time2", "space2", "long_x",
         "basis8", "basis10", "basis12_pre", "basis16", "basis12_filter8", "basis12_filter96",
         "time8_basis24", "space8_basis24", "basis12_two", "two_w09", "two_w10", "two_w11",
         "two_w13", "two_w15", "three_w10", "three_pre_w10", "two_filter8_w10"]
lines += ["", "## Screening evidence", "",
          f"{len(configs)} configurations were screened sequentially at mass "
          "0.1 with 12 RHS, using the same physical operator and full FP64 "
          "checks. These were single samples without a timing warmup; the "
          "selection used operator work and was confirmed in large-batch "
          "finalist/final runs. Screening times are not mixed into the final "
          "timing table. Representative configurations:", "",
          "| Configuration | Aggregate | Vectors/chirality | Filter sweeps | Pre/post | Cycle weight | Coarse dimension | Fine S calls |",
          "|---|---|---:|---:|---:|---:|---:|---:|"]
for name in names:
    r = configs[name]
    lines.append(f"| {name} | {r['aggregate']} | {r['vectors']} | {r['setup_sweeps']} | {r['pre_smooth']}/{r['post_smooth']} | {r['omega']} | {r['coarse_dim']} | {r['fine_schur_calls']} |")
lines += ["", "Removing one smoother alone increased iterations enough to "
          "lose the saved Schur work. A 12-vector basis captured the troublesome "
          "subspace better; 16 vectors did not reduce these counts further. "
          "Shorter time aggregates helped the six-vector basis but did not beat "
          "the twelve-vector 4⁴ configuration. Increasing filtering beyond "
          "eight sweeps did not improve the selected large-batch convergence. "
          "Raising the two-step smoothing weight helped; 1.5 was substantially "
          "worse. These are measured choices for this operator, not universal "
          "multigrid parameter recommendations.", "",
          "## Validation, reproduction, and limits", "",
          f"All {len(rows)} final measured solves passed every-RHS full CPU "
          "FP64 acceptance. Worst final relative residual: "
          f"{max(float(r['max_true_residual']) for r in rows):.3g}; maximum "
          "known-solution relative error: "
          f"{max(float(r['solution_error']) for r in rows):.3g}. The tuned "
          "24-solve small self-test passes, including zero RHS locking. "
          "Both masses' actual 384-dimensional tuned and anisotropic hierarchies "
          "pass P†P, gamma5 Hermiticity, Galerkin consistency, coarse inverse "
          "consistency, CPU/Metal restriction/prolongation, and complete cycle "
          "agreement. Pre-only/post-only and larger-basis paths also have saved "
          "small/full FP64-checked solves.", "",
          "```sh", "make mg-tune-screen", "make mg-tuned-bench", "```", "",
          "The executable's defaults preserve the original settings. To use "
          "the selected cycle directly, add `--vectors 12 --omega 1`; see "
          "[MULTIGRID.md](MULTIGRID.md) for all options and algorithms. "
          "The factory still uses synthetic links and positive masses. The "
          "basis is filtered aggregation rather than a production adaptive "
          "eigensolver; the coarse matrix uses a CPU FP64 direct inverse, "
          "limited to 512 unknowns, with FP32 CPU/GPU applications. Larger "
          "lattices, physical ensembles, and near-critical masses need a "
          "sparse or iterative coarse solve and separate tuning. No "
          "SuperbBLAS or Chroma source is modified.", "",
          "Raw results: [comparison CSV](results/mg-tuned/comparison.csv), "
          "[full FP64 traces](results/mg-tuned/comparison.trace.csv), "
          "[final run log](results/mg-tuned-bench.log), "
          "[screening](results/mg-screen/). The "
          "[metadata note](results/MG_METADATA_NOTE.md) documents correction "
          "of a CSV pre-smoothing metadata field from the recorded run "
          "arguments; numerical measurements were unaffected.", ""]
Path("MULTIGRID_TUNING.md").write_text("\n".join(lines))
print("Wrote MULTIGRID_TUNING.md:", len(rows), "measured solves,", len(configs), "screened configurations")
