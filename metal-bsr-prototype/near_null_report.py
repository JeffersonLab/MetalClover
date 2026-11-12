#!/usr/bin/env python3
"""Validate the paired real-field results and summarize basis/setup tradeoffs."""
import csv
import json
import math
from pathlib import Path
import re
import statistics
import sys


directory = Path(sys.argv[1] if len(sys.argv) > 1 else "results/near-null-svd24")
settings = json.loads((directory / "settings.json").read_text())
expected = set(range(len(settings["rhs_tags"])))
methods = {}
for name, prefix, suffix, method in [
    ("near-null-svd", "near", "mg", "metal-mg"),
    ("filtered-mg", "filtered", "mg", "metal-mg"),
    ("metal-bicgstab", "bicg", "bicg", "metal-bicgstab"),
]:
    with (directory / f"{prefix}.{suffix}.csv").open() as source:
        rows = list(csv.DictReader(source))
    with (directory / f"{prefix}.native.csv").open() as source:
        native = list(csv.DictReader(source))
    if (len(rows) != len(expected) or len(native) != len(expected)
            or {int(r["rhs"]) for r in rows} != expected
            or {int(r["rhs"]) for r in native} != expected):
        raise RuntimeError(f"Missing or duplicated RHS for {name}")
    for row in rows + native:
        residual = float(row["full_relative"])
        if (row["method"] != method or not math.isfinite(residual)
                or not 0 <= residual <= 1e-12):
            raise RuntimeError(f"Physical acceptance failed for {name}: {row}")
    logname = "near-null" if prefix == "near" else prefix
    log = (directory / f"{logname}.log").read_text()
    if "PASS all full FP64 residuals" not in log:
        raise RuntimeError(f"Unfinished solver run: {name}")
    operator = float(re.search(r"original_Chroma_operator_match=([\d.eE+-]+)", log)[1])
    if not math.isfinite(operator) or operator > 1e-13:
        raise RuntimeError(f"Original operator mismatch: {name}")
    setup, setup_work = 0.0, 0
    if suffix == "mg":
        if log.count("MG_READY ") != 1:
            raise RuntimeError(f"Expected exactly one reused hierarchy: {name}")
        ready = re.search(r"MG_READY hierarchy_seconds=([\d.eE+-]+) "
                          r"hierarchy_fine_S=(\d+) GPU_GiB=([\d.eE+-]+)", log)
        setup, setup_work, gpu = float(ready[1]), int(ready[2]), float(ready[3])
    else:
        gpu = float(re.search(r"BICG_READY GPU_GiB=([\d.eE+-]+)", log)[1])
    mean = statistics.fmean(float(r["seconds"]) for r in rows)
    methods[name] = {
        "seconds": [float(r["seconds"]) for r in rows],
        "mean_seconds_per_rhs": mean,
        "hierarchy_seconds": setup,
        "hierarchy_fine_schur_applications": setup_work,
        "measured_seconds_per_rhs_including_setup": mean + setup / len(rows),
        "fine_schur_applications": [int(r["fine_schur_applications"]) for r in rows],
        "measured_fine_schur_per_rhs_including_setup":
            (sum(int(r["fine_schur_applications"]) for r in rows) + setup_work) / len(rows),
        "krylov_iterations": [int(r["krylov"]) for r in rows],
        "defect_corrections": [int(r["outer"]) for r in rows],
        "CPU_FP64_residual_seconds": [float(r["residual_seconds"]) for r in rows],
        "maximum_inline_residual": max(float(r["full_relative"]) for r in rows),
        "maximum_original_Chroma_residual": max(float(r["full_relative"]) for r in native),
        "original_operator_probe": operator,
        "GPU_allocation_at_ready_GiB": gpu,
        "peak_RSS_GiB": float(re.search(r"peak_RSS_GiB=([\d.eE+-]+)", log)[1]),
        "projected_seconds_per_rhs_including_setup": {
            str(n): mean + setup / n for n in (96, 192)
        },
    }

near, old, bi = (methods[k] for k in ("near-null-svd", "filtered-mg", "metal-bicgstab"))
new_mean = near["mean_seconds_per_rhs"]
old_mean = old["mean_seconds_per_rhs"]
bi_mean = bi["mean_seconds_per_rhs"]
null_rows = list(csv.DictReader((directory / "near.near-null.csv").open()))
if len(null_rows) != settings["near_null"]["candidates"]:
    raise RuntimeError("Missing near-null candidate diagnostics")
svd_line = next(line for line in (directory / "near-null.log").read_text().splitlines()
                if line.startswith("NEAR_NULL_SVD "))
svd = {key: float(value) for key, value in
       re.findall(r"(\w+)=([\d.eE+-]+)", svd_line)}
summary = {
    "status": "complete",
    "measured_rhs": len(expected),
    "settings": settings,
    "methods": methods,
    "basis_svd": svd,
    "near_null_candidate_diagnostics": {
        "count": len(null_rows),
        "at_iteration_cap": sum(int(r["iterations"]) == settings["near_null"]["maximum_iterations"]
                                for r in null_rows),
        "minimum_residual_reduction": min(float(r["relative_Sv_to_Sy"]) for r in null_rows),
        "maximum_residual_reduction": max(float(r["relative_Sv_to_Sy"]) for r in null_rows),
    },
    "near_null_vs_bicg_warm_speedup": bi_mean / new_mean,
    "near_null_vs_bicg_speedup_including_measured_setup":
        bi_mean / near["measured_seconds_per_rhs_including_setup"],
    "near_null_vs_filtered_warm_speedup": old_mean / new_mean,
    "fine_work_reduction_vs_bicg":
        1 - sum(near["fine_schur_applications"]) / sum(bi["fine_schur_applications"]),
    "fine_work_reduction_vs_filtered":
        1 - sum(near["fine_schur_applications"]) / sum(old["fine_schur_applications"]),
    "rhs_break_even_vs_filtered":
        max(0, (near["hierarchy_seconds"] - old["hierarchy_seconds"]) / (old_mean - new_mean))
        if old_mean > new_mean else None,
}
(directory / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
for name, row in methods.items():
    print(name, "mean_s=", row["mean_seconds_per_rhs"],
          "setup_s=", row["hierarchy_seconds"],
          "fine_S=", row["fine_schur_applications"],
          "Chroma_max=", row["maximum_original_Chroma_residual"])
print("Near-null / BiCG warm speedup:", summary["near_null_vs_bicg_warm_speedup"])
print("Near-null / filtered warm speedup:", summary["near_null_vs_filtered_warm_speedup"])
