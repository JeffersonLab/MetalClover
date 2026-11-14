# Distillation RHS with Eloy's spin-major caller

One color eigenvector supplies **four physical spinor RHS**, one for each
source spin, per source time and source phasing. The production perambulator
measurement retains Eloy's original ordering: loop over source spin, then
batch eigenvectors within that spin. Its original contractions and storage
organization are restored. The Metal backend keeps its genuine GPU array
implementation and reuses one gauge/clover allocation and MG hierarchy.

| Color eigenvectors (`num_vecs`) | Source spins | Physical RHS |
|---:|---:|---:|
| 1 | 4 | 4 |
| 24 | 4 | 96 |
| 48 | 4 | 192 |
| 5,184 | 4 | 20,736 |

`max_rhs` is the maximum number of physical spinor RHS submitted to qprop
at once. Because each call uses a fixed source spin, it also counts the
number of eigenvectors in that call. **One is valid**, and no rounding to a
multiple of four is applied. A positive cap is respected even when smaller
than four. Nonpositive values retain the helper's existing unlimited behavior.

For example, with `num_vecs=3`, `max_rhs=2`, each spin submits batches of two
and one. Repeating for four spins gives twelve inversions in eight calls.
With `num_vecs=1`, `max_rhs=1`, the four spins make four one-RHS calls.
For 5,184 vectors and `max_rhs=16`, each spin makes 324 calls, giving 1,296
calls and 20,736 physical inversions. `num_vecs` has no small production cap.

`RHSBatch` is an independent maximum resident GPU width (default 8). For
large jobs, `max_rhs=16` and `<RHSBatch>16</RHSBatch>` use the largest width
measured in the existing throughput comparison. A memory-budget estimate
may lower the actual GPU width. A caller cap of one naturally reaches the
scalar GPU path regardless of a larger device cap. The same solver object
and hierarchy survive all four spin loops, source times and source phases.

The temporary full-lattice source/solution pools remain bounded by
`min(num_vecs,max_rhs)`. The localized `Ns` spin-dimension fix is retained
on the Chroma branch. Removing spin coalescing does not remove backend GPU
batching, per-column convergence, nominal work reporting, or native CPU FP64
residual acceptance.

Eloy's measurement contracts and stores one source-spin channel and sink
phasing at a time. It does not retain the overhaul's simultaneous four-spin
projected accumulators. When the eigenvector basis is truncated, sink
projection can substantially reduce stored data. The returned solution tensor,
eigenvectors and projected output still consume memory independently of
`max_rhs`; that parameter bounds inversion temporaries, not every allocation.
The existing callback interface can stream unsmeared output in a separate
writer without restructuring the production perambulator measurement.

## Native throughput harness

`t_metal_clover_batch` takes **BASE_NUM_VECS** and now follows the same
spin-outer ordering as the production measurement. It uses repeatable
complex Z4 color vectors at source time zero. Each vector is embedded into
spins 0–3; logical field index is `4 * vector + spin`, while submission order
is spin first and vector second. CSV output distinguishes `num_vecs`,
`source_spins`, physical `rhs`, seconds/RHS and seconds/colorvec.

Caller caps 1, 8 and 16 are tested, followed by twice as many vectors at cap
16. Actual device width is bounded by the number of vectors in a spin loop:
for one vector all modes still submit one RHS at a time. The base is limited
to 1–24 vectors because this reference-comparison harness retains its source
and solution fields; its final mode reaches 48 vectors / 192 spinor RHS.
This test-only bound does not apply to the production measurement.

```sh
# One vector / four inversions, followed by two vectors / eight inversions:
mainprogs/tests/t_metal_clover_batch gauge.lime action.xml /tmp/eloy1 1

# 24 vectors / 96 inversions, followed by 48 vectors / 192 inversions:
mainprogs/tests/t_metal_clover_batch gauge.lime action.xml /tmp/eloy24 24
```

## Evidence history

The [original throughput report](BATCH_RHS.md) used independent full-field
random spinors. Its 96/192 labels remain physical RHS counts. It demonstrates
the retained Metal backend's batching, rather than a timing of the complete
perambulator job.

The [earlier four-spin-coalescing tests](results/distillation-spins12/summary.json)
record the former implementation at Chroma commit `bc151f436`; they are
historical evidence and do not describe current caller ordering. The restored
spin-major tests record `max_rhs=1`, non-multiple-of-four batches, full native
FP64 residuals and selected direct CPU contraction comparisons separately.


The restored real-gauge regression passed all 64 checks at source time zero.
For both inverters, `max_rhs=1` produced widths `1,1,1,1` for one vector;
`max_rhs=2` produced `2,1` for each spin of three vectors. The new spin-source
cases had maximum full FP64 residual 2.0404e-14 and maximum direct CPU
contraction difference 8.882e-16. The spin-major harness smoke run passed
full-solution comparisons with maximum residual 2.4445e-14 and BiCGStab/MG
difference 2.5393e-13. See [restoration evidence](results/eloy-spin-batching12/summary.json).
