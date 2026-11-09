#!/usr/bin/env python3
"""Validate, then compare original/tuned settings sequentially on one GPU."""
import csv
import os
from pathlib import Path
import subprocess

env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1")
root = Path("results/mg-tuned")
root.mkdir(parents=True, exist_ok=True)


def run(name, args):
    cmd = ["./build/metal_multigrid"] + args
    print("Running", name, " ".join(args), flush=True)
    with (root / (name + ".log")).open("w") as log:
        subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)


run("self-test", ["--self-test", "--vectors", "12", "--omega", "1",
                  "--csv", str(root / "self-test.csv")])
run("hierarchy", ["--validate-only", "--rhs", "3", "--vectors", "12", "--omega", "1"])
run("anisotropic-hierarchy", ["--validate-only", "--rhs", "3", "--aggregate-shape", "4,4,4,2"])

all_rows, all_traces = [], []
cases = [(m, rhs) for m in (".1", ".5") for rhs in ("96", "192")]
for case, (mass, rhs) in enumerate(cases):
    configs = ["original", "tuned"] if case % 2 == 0 else ["tuned", "original"]
    for config in configs:
        name = config + "-m" + mass + "-r" + rhs
        path = root / (name + ".csv")
        args = ["--mass", mass, "--rhs", rhs, "--samples", "3", "--csw", "1",
                "--csv", str(path)]
        if config == "original":
            args += ["--vectors", "6", "--omega", ".65", "--methods",
                     "cpu_bicgstab,metal_bicgstab,cpu_fgmres_mg,metal_fgmres_mg"]
        else:
            args += ["--vectors", "12", "--omega", "1", "--methods",
                     "cpu_fgmres_smooth,metal_fgmres_smooth,cpu_fgmres_mg,metal_fgmres_mg"]
        run(name, args)
        for filename, dest in [(path, all_rows), (Path(str(path) + ".trace.csv"), all_traces)]:
            with filename.open() as f:
                rows = list(csv.DictReader(f))
            for row in rows:
                row["configuration"] = config
                dest.append(row)
        print("Finished", name, flush=True)

for path, rows in [(root / "comparison.csv", all_rows), (root / "comparison.trace.csv", all_traces)]:
    with path.open("w") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
print("Completed", len(all_rows), "measured full FP64-checked solves", flush=True)
