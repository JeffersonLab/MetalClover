#!/usr/bin/env python3
"""Sequential regressions and the 16^4 CPU/Metal benchmark."""
import os
from pathlib import Path
import subprocess

env = dict(os.environ, VECLIB_MAXIMUM_THREADS="1")
Path("results").mkdir(exist_ok=True)
jobs = [
    ("sparse-clover-regression", "metal_clover", ["--self-test"]),
    ("sparse-dense-mg-regression", "metal_multigrid", ["--self-test", "--vectors", "12", "--omega", "1"]),
    ("sparse-self-test-final", "metal_sparse_multigrid", ["--self-test"]),
    ("sparse-hierarchy8-final", "metal_sparse_multigrid", ["--sizes", "8", "--rhs", "3", "--tile-rhs", "3", "--validate-only"]),
    ("sparse-hierarchy16-final", "metal_sparse_multigrid", ["--sizes", "16", "--rhs", "3", "--tile-rhs", "3", "--validate-only"]),
    ("sparse-mg16", "metal_sparse_multigrid", ["--sizes", "16", "--rhs", "96,192", "--tile-rhs", "96", "--mass", ".1,.5", "--coarse-iters", "4", "--samples", "3", "--csw", "1"]),
]
for name, binary, args in jobs:
    cmd = ["./build/" + binary] + args + ["--csv", "results/" + name + ".csv"]
    print("Running:", " ".join(cmd), flush=True)
    with Path("results", name + ".log").open("w") as log:
        subprocess.run(cmd, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)
    print("Completed:", name, flush=True)
print("All sparse study runs passed.", flush=True)
