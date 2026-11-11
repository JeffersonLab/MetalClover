# MetalClover

Standalone Apple GPU experiments for anisotropic clover-Wilson inversion.

The source and reproducible numerical results are in [metal-bsr-prototype](metal-bsr-prototype/).

- [Real 32³×256 low-memory LDL clover, tuned MG timings, and native FP64 checks](metal-bsr-prototype/LOW_MEMORY.md)
- [Real 24³×256 MG tuning and original Chroma FP64 validation](metal-bsr-prototype/MG_TUNING_REAL.md)
- [Native Chroma adapter, action conventions, and initial real-field timings](metal-bsr-prototype/PERFORMANCE_REAL_FIELDS.md)
- [Standalone prototype overview](metal-bsr-prototype/README.md)

Large gauge fields, generated binaries, and saved solution fields are excluded from Git;
text results, input/output hashes, source, and reproduction instructions are versioned.

Run `make format` from the repo root or `metal-bsr-prototype/` to format the
`.hpp`, `.cc`, `.mm`, and `.metal` sources using the checked-in
[`.clang-format`](.clang-format). This requires `clang-format` on `PATH`;
override the executable with `make format CLANG_FORMAT=clang-format-23` if needed.
