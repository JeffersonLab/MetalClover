# Build, install, and checks

Run these commands from the repository root unless a command explicitly changes
directories. The root Makefile defaults to the library; the standalone Makefile
defaults to its `all` program set. Existing experiment-specific instructions
remain in [the prototype README](../metal-bsr-prototype/README.md) and the
individual study documents. This guide collects the supported entry points.

## Requirements

- Apple Silicon macOS with a usable Metal GPU for execution.
- Apple Clang and macOS SDK through Xcode tools (`xcrun clang++`). Offline shader
  checks additionally require `xcrun -sdk macosx metal`, available with the
  appropriate Xcode Metal toolchain.
- `make`, Python 3 (shader embedding and report scripts).
- System Metal, Foundation and Accelerate frameworks; no external BLAS package.
- Optional `clang-format` 23.1.5 to match repository/CI formatting.
- Optional Chroma/QDP++ installation and matching source for native tools.

The build uses C++17/Objective-C++, ARC, blocks, `-O3`, `-mcpu=apple-m1`, and
Accelerate's new LAPACK interface. The library target sets a minimum macOS
deployment version of 15.0 by default (`METALCLOVER_MACOS_MIN`). CI records
build coverage with Xcode 16.4 on macOS 15; that is build evidence, not a
guarantee for every SDK/GPU combination. No CUDA or SuperbBLAS dependency is
needed for the standalone library. Real-field runs require separately supplied
gauge/action data; large binary inputs are not versioned.

Check your tools:

```sh
xcrun --find clang++
xcrun --show-sdk-path
python3 --version
```

## Public library

```sh
make library
make install PREFIX="$HOME/.local/metalclover"
```

The build creates `metal-bsr-prototype/build/api/libmetalclover.dylib` and embeds
shader source into a generated header. Installation writes:

```text
PREFIX/include/metalclover/metalclover.h
PREFIX/lib/libmetalclover.dylib
```

Choose a writable prefix to avoid administrative installation. `DESTDIR` is
also supported for staging. Use `make library METALCLOVER_MACOS_MIN=15.0` to
set the deployment target explicitly. The install name is
`@rpath/libmetalclover.dylib`; clients need an appropriate runtime search path.

For a C++ client with `main` in `client.cc`:

```sh
METALCLOVER_PREFIX="$HOME/.local/metalclover"
xcrun clang++ -std=c++17 -O2 client.cc \
  -I"$METALCLOVER_PREFIX/include" -L"$METALCLOVER_PREFIX/lib" \
  -lmetalclover -Wl,-rpath,"$METALCLOVER_PREFIX/lib" -o client
./client
```

A C client can compile/link with `xcrun clang -std=c11` and the same include,
library and rpath flags. Clients include only the public header; Objective-C++
flags and framework compilation belong to the dylib build. See [API.md](API.md)
for C++ namespacing, data preparation and a compilable function example.

Build the self-contained multi-RHS example against the local library with
`make api-example` from the repository root. `make api-example-run` also runs
its five-source BiCGStab consistency check on the GPU; no gauge files are needed.
See [the runnable example instructions](API.md#c-multi-rhs-example) for MG and
RHS/tile arguments. The executable is `metal-bsr-prototype/build/api/multi_rhs`.

## Standalone programs

```sh
make -C metal-bsr-prototype all
```

This builds nine default programs: `metal_bsr`, `metal_solve`, `metal_wilson`,
`metal_evenodd`, `metal_clover`, `metal_multigrid`, `metal_sparse_multigrid`,
`metal_solver_compare`, and `metal_three_level`, under
`metal-bsr-prototype/build/`. Other tools have explicit targets, for example:

```sh
make -C metal-bsr-prototype build/low_memory_solve build/compact_probe
make -C metal-bsr-prototype -f Makefile.production build/production_solve
```

The library embeds its shaders; standalone programs read shader files relative
to their working directory. Run standalone commands from `metal-bsr-prototype/`
unless their documentation states otherwise. See the respective study for
arguments and input formats; `all` does not build every optional native tool.

## Checks

Formatting and offline shader compilation do not require a numerical GPU run:

```sh
make format-check CLANG_FORMAT=clang-format
make -C metal-bsr-prototype metal-check
```

Selected numerical checks (require GPU access) are:

```sh
make -C metal-bsr-prototype solve-test
make -C metal-bsr-prototype clover-test
make -C metal-bsr-prototype mg-test
make -C metal-bsr-prototype sparse-mg-test
make -C metal-bsr-prototype three-test
```

These exercise different solver stages. The generic `test` target checks the
BSR prototype only; it is not an umbrella for all solvers. Several targets set
`VECLIB_MAXIMUM_THREADS=1` to avoid BLAS oversubscription. For repeatable timing,
follow the experiment-specific warmup, tiling and validation procedure instead
of treating a self-test time as a benchmark.

The public API contract test requires a small real aligned `CHRMG003` bridge:

```sh
make -C metal-bsr-prototype api-contract-test BRIDGE=/path/to/aligned.bridge
```

Supply the bridge file expected by `api_contract_test.cc`; generated native
input fields and full residual checks are separate from hosted build checks.
This Make target runs the default BiCGStab contract test. To exercise MG with
the same built test program, run from `metal-bsr-prototype/`:

```sh
VECLIB_MAXIMUM_THREADS=1 ./build/api/api_contract_test /path/to/aligned.bridge mg
```

## Optional native Chroma tools and integration

To build native exporters/checkers, point to an existing Chroma installation
and a source snapshot matching that installation:

```sh
make -C metal-bsr-prototype -f Makefile.production production \
  CHROMA_CONFIG=/path/to/chroma/bin/chroma-config \
  CHROMA_SOURCE=/path/to/matching/chroma/source
make -C metal-bsr-prototype -f Makefile.production build/low_memory_verify \
  CHROMA_CONFIG=/path/to/chroma/bin/chroma-config \
  CHROMA_SOURCE=/path/to/matching/chroma/source
```

The production Makefile uses the installed Chroma compiler/flags for native
tools, Apple Clang for Metal programs, and `BASE_PRECISION=64` for native
reference parsing/checking. It does not modify Chroma or SuperbBLAS trees.
Use [the adapter guide](../metal-bsr-prototype/CHROMA_ADAPTER.md) and
[low-memory guide](../metal-bsr-prototype/LOW_MEMORY.md) for bridge generation,
action XML and full verification workflows.

For the separate Chroma branch that implements MetalClover support, install
the library first, then add `--with-metal=/absolute/install/prefix` to Chroma's
existing configure command, retaining its normal QDP++/other dependency flags.
That integration provides `MetalBiCGClover` and `MetalMGClover` and requires
the batch API extension. An arbitrary upstream Chroma checkout may not contain
these adapters or tests. Build native regression programs in the Chroma build
tree with:

```sh
make -C mainprogs/tests t_metal_clover t_metal_clover_batch
```

See [BATCH_RHS.md](../metal-bsr-prototype/BATCH_RHS.md) for invocation and
`<RHSBatch>` controls. The external adapter's complete XML reference lives
in that Chroma checkout's `docs/metal-clover.md`.

## Common build/run issues

| Symptom | Check |
| --- | --- |
| `xcrun` cannot find a tool/SDK | Select an installed Xcode/toolchain with the needed macOS/Metal tools. |
| `dyld` cannot load the library | Link with the installed prefix's `lib` directory in rpath. |
| Standalone shader file not found | Run from `metal-bsr-prototype/`. |
| Invalid MG aggregation geometry | Each fine block must divide its extent and produce an even first coarse extent. |
| Unexpected memory pressure | Lower RHS cap/restart or caller group size; use the documented low-memory path. |
| Native header/link mismatch | Use Chroma sources matching `chroma-config` and the installed library. |
| FP64 convergence failure | Check action, spinor layout/basis, boundary handling and iteration budget; inspect `last_error`. |

Rebuild the library after editing its implementation/header or shader sources.
Keep build products and large generated fields outside the documentation tree.
