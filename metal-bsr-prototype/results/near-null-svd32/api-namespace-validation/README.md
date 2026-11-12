# Post-benchmark interface and namespace validation

The final driver uses `MetalCloverEnv`, the explicit Metal 3.1 shader standard,
and the structured end-of-run work report. It was tested on the real 12³×96
field with two deterministic RHS. Both solutions are bitwise identical to the
pre-refactor control; `native.csv` records independent original-Chroma FP64
full residuals below 1e-12. These are correctness tests performed alongside a
CPU build, not performance comparisons.

The GCC 14 client links the Apple Clang dylib through the installed C ABI.
Its contract check covers full and odd solves, full reconstruction, scalar
versus array solutions, empty arrays, exact initial guesses, invalid iteration
budgets and tolerances, aliased outputs, invalid ABI and invalid local LDL
pivots. `installed-api.json` records the public symbols and installed hashes.

The immutable large benchmark predates this API refactor. Recompiling its
saved source snapshot produced a byte-identical executable; the original
and reconstructed binary hashes are in `benchmark-reproduction.json`.

The six modified native Chroma bridge/checker C++ programs also pass syntax
checks against the matching installed FP64 Chroma source headers; see
`native-cpp-results.json`. All 19 standalone Objective-C++ programs pass
syntax checks. Formatting covers the full `.hpp`, `.cc`, `.mm`, `.metal` and
public `.h` set.

Versioned text logs/CSV use LF line endings and omit trailing whitespace;
measurements are unchanged. The source patch is gzip-compressed without
changing its bytes, with its uncompressed hash in `benchmark-patch.json`.
