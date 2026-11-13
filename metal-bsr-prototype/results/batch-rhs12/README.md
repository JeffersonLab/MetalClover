# Native Chroma batched RHS evidence

[BATCH_RHS.md](../../BATCH_RHS.md) explains the implementation, controls,
measurement scope and results. These are local, single-rank tests on the real
12³×96 LIME field at mass −0.0743, using the original native Chroma FP64
operators and a 1e-12 full-system relative residual tolerance.

- `action.xml` is the exact action metadata used; `inputs.json` hashes it and
  the local gauge field. Neither the gauge field nor binary reference fields
  are stored in Git.
- `run.csv`, `summary.json`, `timing.log`, `status.json` and `resources.csv`
  record the completed 96/192-RHS comparison, actual device widths, work
  counts, one reused hierarchy, acceptance and resource samples.
- `validation/` contains the preceding native Chroma regression, C ABI
  tests, GCC-client check, syntax/format checks, fresh configure checks,
  build logs and local reproduction controllers. The benchmark smoke test
  uses only two/four RHS; its timings are not throughput results.
- `validation/benchmark-source-binary.json` was captured before timing and
  verified unchanged afterward. `validation/sources.json` records the wider
  source inventory after timing, base commits and the exact user style hash.
- `validation/timing-build-log-copies.json` records original and normalized
  hashes for the final timing/build evidence. Copies use LF and trim trailing
  whitespace; original local logs remain at their recorded scratch paths.
- `validation/reference-artifact.json` records the local 96-source scalar
  reference file's size and SHA-256. Complete solutions for those 96 sources
  were compared for every width and both methods. All 192 sources received
  independent full native residual checks.

Hierarchy construction is excluded from solve time and recorded separately.
All native qprop reduction, reconstruction and FP64 acceptance work is inside
solve time. Additional post-timing verification, source generation and
reference-file I/O are outside. Measurements are single runs. Rates count
nominal dispatched GPU operator work, including inactive lanes, divided by
complete qprop wall time; they are not hardware-counter rates.

The matching Chroma integration is
[`b3a380ea6`](https://github.com/JeffersonLab/chroma/commit/b3a380ea62c7c9002781a5bc1e9a557e4d3f3521)
on `rge/devel-datalib-reorg`. The MetalClover implementation is in the commit
containing this evidence directory. `validation/integration.json` records the
branch pairing and final formatting checks.
