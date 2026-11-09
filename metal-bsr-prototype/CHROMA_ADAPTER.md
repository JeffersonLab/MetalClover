# Standalone Chroma adapter

`chroma_export` uses the installed Chroma and QDP++ libraries to read a NERSC
archive. It writes the forward links in x-fastest order and two FP64 reference
actions. The exporter constructs the full operator from Chroma's
`QDPCloverTermD` and `QDPWilsonDslashD` and builds the unscaled odd Schur
reference using Chroma's even clover inverse. The standalone executable
reconstructs its own clover blocks and checks both actions against the bridge
before building the hierarchy. The reference probe covers every spin, color,
site, and boundary seam.

The adapter handles single-rank isotropic SU(3) with four spin components,
an even hypercube with extent 4–32, `c_sw=1`, `u0=1`, and fermion phases
`(1,1,1,-1)`. The finest aggregation size must divide the extent and leave
an even first coarse extent. The Chroma mass travels in the bridge, so the
solver takes the lattice extent and mass from that file. The archived gauge
field is loaded at the installed QDP++ build's base precision. Both FP64
operators use the same loaded links promoted to double.

## Build and run

```sh
make build/chroma_mg build/chroma_export \
  CHROMA_CONFIG=/path/to/chroma/bin/chroma-config
```

Set `CHROMA_EXTRA_CPPFLAGS=-I/path/to/matching/chroma/lib` if the installed
headers omit `coarsenable_op.h`. Use the source tree that matches the installed
headers and library. The locally installed Chroma 3.44.0 used
`/Users/admin/Downloads/qcd-macos-share-20261005/chroma.devel/lib` for
that header.

```sh
OMP_NUM_THREADS=1 ./build/chroma_export /path/to/nersc.config build/field.bridge -.375
VECLIB_MAXIMUM_THREADS=1 ./build/chroma_mg \
  --chroma-field build/field.bridge --block 5 --tile-rhs 4 \
  --rhs 96,192 --methods mg3 --csw 1 --threads 8 --check-hierarchy \
  --csv results/chroma-run.csv
```

The example uses a 20⁴ field. A single run builds the two sparse coarse
matrices once and reuses them across both backends and all RHS tiles.
The 96-RHS timing is the first 24 measured tiles and the 192-RHS timing
includes all 48. Supply a fresh CSV filename for an imported-field run.

## RHS and solutions

The benchmark generates deterministic Z4 sources directly in each tile.
For caller-supplied spinors, pass `--chroma-rhs FILE`. The file contains
exactly the largest requested RHS count. Layout is RHS, x-fastest site,
color, spin, then real/imaginary IEEE binary64 in little-endian order.
For each lexicographic site coordinate, obtain the Chroma local index via
`QDP::Layout::linearSiteIndex(coord)`. Read a component from
`field.elem(index).elem(spin).elem(color)` and write its real and imaginary
parts in that order.

Pass `--chroma-solutions PREFIX` with a single MG sample to write
`PREFIX.metal.bin` and/or `PREFIX.cpu.bin` in the same layout. Each tile is
written after all its full residuals pass. Input transposition and output
packing are included in the recorded tile wall time when used. The prototype
keeps only one RHS tile in memory.

```sh
OMP_NUM_THREADS=1 ./build/chroma_export --check-solutions \
  /path/to/nersc.config PREFIX.metal.bin -.375 /path/to/rhs.bin
```

Omit the final RHS path when checking the generated Z4 sources. This mode
loads each returned solution into Chroma and independently checks the full
FP64 residual. The solver itself checks the full CPU FP64 residual for every
RHS after every outer correction, with relative 2-norm target `1e-10`.

The bridge is a local interchange file with magic `CHRMG001`, a little-endian
32-bit extent, an FP64 mass, then the complex FP64 forward links ordered by
site, direction, row, and column. The full and odd-only reference vectors
follow in site, color, spin order. Its convention and version are fixed by
`chroma_bridge.hpp`.
