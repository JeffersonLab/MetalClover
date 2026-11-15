from pathlib import Path
import json
import sys

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
t = int(sys.argv[1])
report = ROOT / f'metal-bsr-prototype/results/perambulator-t{t}-5184'
point = json.loads((report / 'point.complete.json').read_text())
run = json.loads((report / 'point.run.json').read_text())
peram = json.loads((report / 'peram.complete.json').read_text())
pipeline = json.loads((report / 'point.pipeline.json').read_text())
staged = json.loads((report / 'staged-run.json').read_text())
validation = json.loads((report / 'production-validation.json').read_text())
assert validation['status'] == 'pass'


def duration(seconds):
    rounded = round(seconds)
    return f'{rounded // 3600} h {rounded // 60 % 60} min {rounded % 60} s'


text = f'''# Configuration 1000: standard perambulator at source/sink t={t}

**Completed and verified, {pipeline['finished_UTC'][:10]} UTC.** All 20,736 RHS
passed both original full Chroma CPU FP64 residual checks at 1e-12. One MG
hierarchy served all 1,296 width-sixteen array calls in one uninterrupted session.
Native rereading checked every saved point column; standard conversion passed
all sixteen full record readbacks, direct QDP contractions, full gamma5
Hermiticity and the independent DataLib key check.

| Measurement | Result |
| --- | ---: |
| Solve-session wall, including input/setup/full checks/output | {point['elapsed_seconds_across_sessions']:,.3f} s ({duration(point['elapsed_seconds_across_sessions'])}) |
| Native qprop calls, including native full FP64 checks | {point['qprop_seconds']:,.3f} s |
| Additional independent full FP64 checks | {point['verify_seconds']:,.3f} s |
| Point extraction/write/readback | {point['write_seconds']:.3f} s |
| Maximum full relative residual, either check | {point['maximum_full_relative']:.16g} |
| Iterations per physical RHS | {point['iterations_min']}–{point['iterations_max']} |
| Maximum sampled solve process RSS | {run['peak_RSS_GiB']:.3f} GiB |
| Nominal dispatched GPU operator work | {point['nominal_dispatched_operator_flops']:,} FLOP |
| Nominal GPU operator work / total solve wall | {point['nominal_operator_GFLOP_per_second_total_wall']:.2f} GFLOP/s |
| Standard projection plus complete native validation | {peram['native_seconds']:.3f} s |
| Complex-double GEMM time | {peram['gemm_seconds']:.3f} s |
| Maximum sampled projection process RSS | {peram['peak_RSS_GiB']:.3f} GiB |
| Full projected gamma5 relative Frobenius discrepancy | {peram['full_gamma5_relative_Frobenius']:.15g} |
| Maximum absolute gamma5 discrepancy | {peram['full_gamma5_max_absolute']:.15g} |
| Maximum independent QDP contraction error, norm-scaled | {peram['max_QDP_dot_scaled_error']:.15g} |
| Complete local pipeline, including whole-file hashes and reader | {pipeline['elapsed_seconds']:,.3f} s ({duration(pipeline['elapsed_seconds'])}) |

Hardware: Apple M1 Max, 24 GPU cores, 32 GiB, one rank. The nominal GPU
operator numerator includes dispatched inactive/padded lanes and excludes CPU
FP64 and Krylov-vector arithmetic. Its denominator includes the complete
solve-session wall; final checksum rereading/hashing and standard projection
are separate. Process RSS is sampled every five seconds.

The verified standard perambulator is:

```
{staged['standard_perambulator']}
```

It contains sixteen 5,184×5,184 complex-FP64 matrices with source and sink time
{t} in the actual keys and XML, mass label `light`, and the original vector
ordering and phases. Size is {peram['output_bytes']:,} bytes; SHA256:
`{peram['output_SHA256']}`.

The accepted intermediate point file is `{staged['point_prefix']}.npy`,
{point['payload_bytes']:,} bytes, SHA256
`{point['payload_sha256']}`.
It retains the unsmeared solution only at the common sink time after the full
four-dimensional native qprop and independent original Chroma FP64 checks.

The supplied eigenvectors are `{staged['eigenvector_file']}`,
{staged['eigenvector_bytes']:,} bytes, SHA256
`{staged['eigenvector_SHA256']}`. Their header reports 12³,
decay_dir=3, t_slice={t} and 5,184 vectors. Stored FP32 values are promoted to
FP64; the residual checks refer to those supplied promoted sources.

Only `Contractions/t_sources` changes to {t} in `action.xml`. The mass −0.0743,
anisotropies, clover coefficients, spatial stout action and antiperiodic
temporal boundary conditions are identical to the earlier accepted runs.
No numerical source, executable, library or solver setting changed.

`production-validation.json` checks every source-spin collection, all native
array calls, one hierarchy/session, both residual bounds, complete file checks
and frozen source/executable/library hashes. `point.*` and `peram.*` contain
the full journals, metadata, logs, checksums and resource samples. The queue
and pipeline scripts run locally with sleep prevention, independently of the
user's iPhone connection. Large output data remain outside Git.
'''
if t == 2:
    text += '''
The successful time-two pipeline automatically triggered the queued time-three
repeat, including the dependency on standard validation and independent reader
completion. Its separate hierarchy is reused for that entire new collection.
See [the time-three record](../perambulator-t3-5184/README.md).
'''
text += '''
Storage: with explicit user approval, the earlier time-zero and time-one
unsmeared `.npy` intermediates were removed after checking their hashes and
rechecking the retained standard perambulator hashes. Every standard file,
original gauge/eigenvector input, log and residual record is retained; see
[the cleanup record](../perambulator-t3-5184/intermediate-cleanup.json).
'''
cleanup = ROOT / 'metal-bsr-prototype/results/perambulator-t4-5184/intermediate-cleanup.json'
if cleanup.exists():
    saved = json.loads(cleanup.read_text())
    if saved.get('removed_intermediate_files'):
        text += '''
The same verified cleanup procedure subsequently retired the time-two and
 time-three point intermediates to make room for time four. Their complete
standard perambulators, accepted logs, residuals and hashes are retained; see
[that cleanup record](../perambulator-t4-5184/intermediate-cleanup.json).
'''.replace('\n time-three', '\ntime-three')
cleanup5 = ROOT / 'metal-bsr-prototype/results/perambulator-t5-5184/intermediate-cleanup.json'
if cleanup5.exists() and json.loads(cleanup5.read_text()).get('removed_intermediate_files'):
    text += '''
The accepted time-four point intermediate was also retired after complete
evidence capture and point/standard hash rechecks to make room for time five.
Its standard perambulator and all original inputs/logs remain; see
[the cleanup record](../perambulator-t5-5184/intermediate-cleanup.json).
'''
(report / 'README.md').write_text(text)
print(f'Wrote completed time-{t} report.')
