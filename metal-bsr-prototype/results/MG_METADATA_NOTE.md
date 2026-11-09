# Pre-smoothing metadata correction

The first tuning binary accidentally wrote the preconditioner-call counter into
its new `pre_smooth` CSV column because of a local variable name collision.
The numerical algorithm, all timings, residuals, and operator/cycle counters
were unaffected. The driver now uses distinct names.

For the saved tuning data only, `pre_smooth` was restored from the actual
`pre/post=` settings in each corresponding run log. Combined tables inherit
that corrected metadata. No numerical measurements were edited. Original CSVs
before this metadata correction are preserved in
`mg-pre-smooth-metadata-originals/`.
