from pathlib import Path
import datetime
import hashlib
import json
import shutil
import subprocess

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
REPORT = ROOT / 'metal-bsr-prototype/results/perambulator-t5-5184'
BASE = Path('/Users/admin/Downloads/MetalClover-results')
record = {
    'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'authorization': 'Requested continuation of the same queue procedure: retain every standard perambulator and retire only the accepted previous-run point intermediate after complete evidence and hash verification',
    'minimum_free_disk_GiB_before_start': 17.0,
    'free_disk_GiB_before': shutil.disk_usage(BASE).free / 2**30,
    'all_original_inputs_logs_checksums_and_standard_perambulators_retained': True,
    'removed_intermediate_files': [],
    'retained_standard_perambulators': [],
}


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


if record['free_disk_GiB_before'] < record['minimum_free_disk_GiB_before_start']:
    point = BASE / 'cfg1000-t4-unsmeared'
    payload = Path(str(point) + '.npy')
    peram = BASE / 'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_4.sdb1000'
    prior_report = ROOT / 'metal-bsr-prototype/results/perambulator-t4-5184'
    accepted = read(str(point) + '.complete.json')
    projected = read(str(peram) + '.complete.json')
    proof = read(prior_report / 'production-validation.json')
    assert accepted['status'] == projected['status'] == 'complete'
    assert proof['status'] == 'pass'
    assert accepted['source_time'] == accepted['sink_time'] == 4
    assert projected['source_time'] == projected['sink_time'] == 4
    assert accepted['completed_rhs'] == 20736 and projected['spin_records'] == 16
    assert accepted['maximum_full_relative'] <= 1e-12
    assert accepted['maximum_qprop_relative'] <= 1e-12
    assert projected['full_gamma5_relative_Frobenius'] < 1e-12
    assert proof['point_SHA256'] == accepted['payload_sha256']
    assert proof['standard_perambulator_SHA256'] == projected['output_SHA256']
    for path, size, expected_sha, kind in [
        (payload, accepted['payload_bytes'], accepted['payload_sha256'], 'intermediate'),
        (peram, projected['output_bytes'], projected['output_SHA256'], 'standard'),
    ]:
        assert path.stat().st_size == size, path
        assert sha(path) == expected_sha, path
        entry = {'path': str(path), 'bytes': size, 'accepted_SHA256': expected_sha,
                 'current_SHA256_matches_accepted_result': True}
        record['removed_intermediate_files' if kind == 'intermediate'
               else 'retained_standard_perambulators'].append(entry)
    assert payload == BASE / 'cfg1000-t4-unsmeared.npy'
    opened = subprocess.run(['/usr/sbin/lsof', '-Fpn', '--', str(payload)],
                            text=True, capture_output=True)
    assert opened.returncode == 1 and not opened.stdout, opened.stdout + opened.stderr
    record['intermediate_files_not_open'] = True
    payload.unlink()
    record['status'] = 'reclaimed_verified_intermediate'
    proof.update(point_intermediate_retired_after_full_validation=True,
                 retirement_record='../perambulator-t5-5184/intermediate-cleanup.json',
                 standard_perambulator_retained_SHA256_rechecked=True)
    (prior_report / 'production-validation.json').write_text(json.dumps(proof, indent=2) + '\n')
    readme = prior_report / 'README.md'
    readme.write_text(readme.read_text() + '''

Storage update for the time-five repeat: the accepted time-four unsmeared
`.npy` intermediate was retired after complete evidence capture, point/standard
SHA256 rechecks and confirmation that it was not open. Its standard perambulator,
original inputs and all logs/residual records are retained. See
[the cleanup record](../perambulator-t5-5184/intermediate-cleanup.json).
''')
else:
    record['status'] = 'no_cleanup_needed'
record['free_disk_GiB_after'] = shutil.disk_usage(BASE).free / 2**30
(REPORT / 'intermediate-cleanup.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
