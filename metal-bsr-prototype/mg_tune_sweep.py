#!/usr/bin/env python3
"""Sequential, reproducible screening; never run GPU benchmarks concurrently."""
import csv
import os
from pathlib import Path
import subprocess
import sys

configs = [
    ("original", 6, 8, "4,4,4,4", 1, 1, .65, 8),
    ("pre_only", 6, 8, "4,4,4,4", 1, 0, .65, 8),
    ("post_only", 6, 8, "4,4,4,4", 0, 1, .65, 8),
    ("filter24", 6, 24, "4,4,4,4", 0, 1, .65, 8),
    ("filter48", 6, 48, "4,4,4,4", 0, 1, .65, 8),
    ("basis12", 12, 24, "4,4,4,4", 0, 1, .65, 8),
    ("basis12_filter48", 12, 48, "4,4,4,4", 0, 1, .65, 8),
    ("time2", 6, 24, "4,4,4,2", 0, 1, .65, 8),
    ("time2_filter48", 6, 48, "4,4,4,2", 0, 1, .65, 8),
    ("space2", 6, 24, "4,4,2,4", 0, 1, .65, 8),
    ("long_x", 6, 24, "8,4,4,2", 0, 1, .65, 8),
    ("long_x_basis12", 12, 24, "8,4,4,2", 0, 1, .65, 8),
    ("basis12_pre", 12, 24, "4,4,4,4", 1, 0, .65, 8),
    ("time2_pre", 6, 24, "4,4,4,2", 1, 0, .65, 8),
    ("time2_w09", 6, 24, "4,4,4,2", 0, 1, .9, 8),
    ("time2_w11", 6, 24, "4,4,4,2", 0, 1, 1.1, 8),
    ("time2_restart12", 6, 24, "4,4,4,2", 0, 1, .9, 12),
    ("basis12_w09", 12, 24, "4,4,4,4", 0, 1, .9, 8),
]
if "--stage2" in sys.argv:
    configs = [
        ("basis8", 8, 24, "4,4,4,4", 1, 0, .65, 8),
        ("basis10", 10, 24, "4,4,4,4", 1, 0, .65, 8),
        ("basis16", 16, 24, "4,4,4,4", 1, 0, .65, 8),
        ("basis12_filter8", 12, 8, "4,4,4,4", 1, 0, .65, 8),
        ("basis12_filter16", 12, 16, "4,4,4,4", 1, 0, .65, 8),
        ("basis12_filter96", 12, 96, "4,4,4,4", 1, 0, .65, 8),
        ("basis12_two", 12, 24, "4,4,4,4", 1, 1, .65, 8),
        ("basis12_two_pre", 12, 24, "4,4,4,4", 2, 0, .65, 8),
        ("basis16_two", 16, 24, "4,4,4,4", 1, 1, .65, 8),
        ("time8_basis24", 24, 24, "4,4,4,8", 1, 0, .65, 8),
        ("space8_basis24", 24, 24, "8,4,4,4", 1, 0, .65, 8),
        ("basis12_restart12", 12, 24, "4,4,4,4", 1, 0, .65, 12),
    ]
if "--stage3" in sys.argv:
    configs = [
        ("two_w09", 12, 24, "4,4,4,4", 1, 1, .9, 8),
        ("two_w10", 12, 24, "4,4,4,4", 1, 1, 1., 8),
        ("two_w11", 12, 24, "4,4,4,4", 1, 1, 1.1, 8),
        ("two_w13", 12, 24, "4,4,4,4", 1, 1, 1.3, 8),
        ("two_w15", 12, 24, "4,4,4,4", 1, 1, 1.5, 8),
        ("basis6_two_w10", 6, 24, "4,4,4,4", 1, 1, 1., 8),
        ("basis8_two_w10", 8, 24, "4,4,4,4", 1, 1, 1., 8),
        ("time2_two_w10", 6, 24, "4,4,4,2", 1, 1, 1., 8),
        ("time2_two_w13", 6, 24, "4,4,4,2", 1, 1, 1.3, 8),
        ("three_w10", 12, 24, "4,4,4,4", 2, 1, 1., 8),
        ("three_pre_w10", 12, 24, "4,4,4,4", 3, 0, 1., 8),
        ("two_filter8_w10", 12, 8, "4,4,4,4", 1, 1, 1., 8),
    ]

env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1")
out = Path("results/mg-screen")
out.mkdir(parents=True, exist_ok=True)
all_rows = []
for name, nv, sweeps, shape, pre, post, omega, restart in configs:
    path = out / (name + ".csv")
    cmd = ["./build/metal_multigrid", "--rhs", "12", "--mass", ".1",
           "--vectors", str(nv), "--setup-sweeps", str(sweeps),
           "--aggregate-shape", shape, "--pre-smooth", str(pre),
           "--post-smooth", str(post), "--omega", str(omega),
           "--restart", str(restart), "--samples", "1", "--no-warmup",
           "--methods", "cpu_fgmres_mg,metal_fgmres_mg", "--csv", str(path)]
    with path.with_suffix(".log").open("w") as log:
        subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    with path.open() as f:
        rows = list(csv.DictReader(f))
    for row in rows:
        row["configuration"] = name
        all_rows.append(row)
        print(name, row["backend"], "inner", row["inner_max_sum"],
              "fineS", row["fine_schur_calls"], "ms", row["total_ms"], flush=True)
stage = "3" if "--stage3" in sys.argv else "2" if "--stage2" in sys.argv else ""
with (out / ("screen" + stage + ".csv")).open("w") as f:
    writer = csv.DictWriter(f, fieldnames=list(all_rows[0]))
    writer.writeheader()
    writer.writerows(all_rows)
