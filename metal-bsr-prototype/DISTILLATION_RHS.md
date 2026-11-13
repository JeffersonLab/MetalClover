# Four-spin distillation RHS

A distillation eigenvector is a spatial color vector. For each source time
and source phasing, it produces **four physical spinor RHS**, one for each
source spin. The same color vector is embedded separately into spins 0–3.
Sink spin is an output index and does not multiply the number of inversions.

| Color eigenvectors (`num_vecs`) | Source spins | Physical RHS |
|---:|---:|---:|
| 1 | 4 | 4 |
| 2 | 4 | 8 |
| 24 | 4 | 96 |
| 48 | 4 | 192 |
| 96 | 4 | 384 |

The Chroma distillation measurement submits complete eigenvector groups to
the array qprop interface. The flattened index is
`rhs = 4 * vector_index + source_spin`. This makes four independent solves
available together even when `num_vecs=1`. The gauge field, clover data and
one MG hierarchy are reused across vectors, spins, source times and phases
served by the same solver object.

`max_rhs` continues to count **physical spinors**, rather than eigenvectors.
The measurement normalizes a positive cap down to a multiple of four, with
a minimum of four, and bounds it by the requested total. Thus caps 1–7
select four RHS; 8–11 select eight. A nonpositive cap selects the full
four-spin collection. `num_vecs` itself is never rounded or increased.
Logs show the requested/effective cap, vectors per batch and physical RHS
per source time and phasing. Total inversion count is
`4 * num_vecs * number_of_source_times * number_of_source_phasings`.
Additional sink phasings reuse those inversions.

`<RHSBatch>16</RHSBatch>` independently caps resident GPU RHS. The generic
Metal library still accepts arbitrary RHS arrays and may split a caller's
four-spin group to respect device memory. Setting this GPU cap to one is
also supported for a scalar performance baseline. At the measurement level,
`max_rhs=16` supplies up to four complete eigenvectors per qprop call.

The measurement retains projected perambulator elements while streaming
bounded full-field tiles. Source vector, source spin and sink spin remain
separate tensor/storage indices. All per-column checks use the original
Chroma CPU FP64 operator at a 1e-12 relative residual target.

The current native throughput harness takes **BASE_NUM_VECS**, not a spinor
RHS count. It uses repeatable time-localized complex Z4 color vectors at
`t_source=176 mod Lt`, embedding each in all four spins. CSV output includes
`num_vecs`, `source_spins`, physical `rhs`, seconds/RHS and seconds/colorvec.
It supports a base of 1–24 vectors and doubles the collection in its last
measurement, keeping the original 192-spinor memory bound.

```sh
# 1 vector / 4 RHS, followed by 2 vectors / 8 RHS:
mainprogs/tests/t_metal_clover_batch gauge.lime action.xml /tmp/spins1 1

# 24 vectors / 96 RHS, followed by 48 vectors / 192 RHS:
mainprogs/tests/t_metal_clover_batch gauge.lime action.xml /tmp/spins24 24
```

The earlier [batch timing report](BATCH_RHS.md) used independent full-field
random spinors. Its 96/192 labels are physical RHS counts; those timings
have not been relabeled as eigenvector or four-spin distillation timings.

## Validation

On the real 12³×96 gauge field at mass −0.0743, the native regression passed
64 checks. New four-spin cases used one vector / four RHS in one GPU tile,
and three vectors / twelve RHS in tiles of eight plus four. Both methods
passed full original Chroma CPU FP64 residual checks at 1e-12; the maximum
for these spin-source cases was 2.2152e-14. Selected projected elements at
the source time agreed with direct CPU inner products within 8.882e-16.

The updated throughput harness also passed its one/two-vector smoke test,
including scalar/batch and complete BiCGStab/MG solution comparisons. It
reported the required factor of four between seconds/colorvec and
seconds/RHS, with maximum full residual 2.1274e-14 and solution difference
4.0540e-13. These small-run timings are not a new 96/192-RHS performance
comparison. This validates real-gauge qprop, source assembly and tiled
contractions; a complete inline perambulator file-I/O job was not run.

[Validation records](results/distillation-spins12/summary.json) contain the
accepted counts, residuals, projections, source hashes and smoke CSV.
