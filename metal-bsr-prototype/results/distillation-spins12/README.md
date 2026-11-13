# Four-spin distillation validation

See [DISTILLATION_RHS.md](../../DISTILLATION_RHS.md) for the source layout,
normalization rules and results. The real 12³×96 LIME field and exact action
metadata are hashed in `inputs.json`; gauge and binary reference fields are
excluded from Git.

`regression.log` records 64 accepted native checks, including 32 new four-spin
cases. One colorvec supplies four physical RHS in one GPU tile; three supply
twelve in tiles of eight and four. Every full spin-source solution is checked
with the original Chroma CPU FP64 operator. Selected source-time perambulator
elements are checked against direct CPU inner products.

`spins1.csv` and `four-spin-smoke.log` record one/two-colorvec harness runs:
four/eight physical RHS, four source spins, scalar and GPU-array qprop calls,
one reused MG hierarchy, full FP64 checks and full-solution comparisons.
Small smoke times are not 96/192-RHS throughput measurements.

`validation-status.json` records final build/run acceptance. The changed
native objects and dependent Chroma routines were rebuilt before relinking.
`sources.json` hashes the final canonical/tested source copies. `artifacts.json`
hashes the actual binaries, dylib and local scalar reference file. Build logs in `compile-logs/`,
reproduction controllers and copy-normalization hashes are included. Copied
logs use LF and trim trailing whitespace. This tests qprop, colorvec source
assembly and the tiled contraction pattern; a complete inline measurement
with perambulator file I/O was not run.

The matching Chroma change is
[`bc151f436`](https://github.com/JeffersonLab/chroma/commit/bc151f436806191452d3aa363ae4ca2972e3ceaf)
on `rge/devel-datalib-reorg`. The MetalClover evidence/documentation update is
in the commit containing this directory. `integration.json` records the
branch pairing and final checks. The MetalClover numerical library is unchanged
from the preceding batch implementation.
