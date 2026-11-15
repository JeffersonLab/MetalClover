from pathlib import Path
import datetime
import hashlib
import json
import shutil
import subprocess

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
REPORT = ROOT / 'metal-bsr-prototype/results/perambulator-t4-5184'
BASE = Path('/Users/admin/Downloads/MetalClover-results')
record = {
    'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'authorization': 'User requested time four, same as before; retain all standard perambulators and reclaim only validated point intermediates as in the explicitly approved time-three procedure',
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
    intermediates = []
    for t in (2, 3):
        point = BASE / f'cfg1000-t{t}-unsmeared'
        peram = BASE / f'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb1000'
        accepted = read(str(point) + '.complete.json')
        projected = read(str(peram) + '.complete.json')
        proof = read(ROOT / f'metal-bsr-prototype/results/perambulator-t{t}-5184/production-validation.json')
        assert accepted['status'] == projected['status'] == 'complete'
        assert proof['status'] == 'pass'
        assert accepted['source_time'] == accepted['sink_time'] == t
        assert projected['source_time'] == projected['sink_time'] == t
        assert accepted['completed_rhs'] == 20736 and projected['spin_records'] == 16
        assert accepted['maximum_full_relative'] <= 1e-12
        assert accepted['maximum_qprop_relative'] <= 1e-12
        assert projected['full_gamma5_relative_Frobenius'] < 1e-12
        assert proof['point_SHA256'] == accepted['payload_sha256']
        assert proof['standard_perambulator_SHA256'] == projected['output_SHA256']
        for path, expected_size, expected_sha, kind in [
            (Path(str(point) + '.npy'), accepted['payload_bytes'], accepted['payload_sha256'], 'intermediate'),
            (peram, projected['output_bytes'], projected['output_SHA256'], 'standard'),
        ]:
            assert path.stat().st_size == expected_size, path
            assert sha(path) == expected_sha, path
            entry = {'path': str(path), 'bytes': expected_size,
                     'accepted_SHA256': expected_sha,
                     'current_SHA256_matches_accepted_result': True}
            if kind == 'intermediate':
                intermediates.append(path)
                record['removed_intermediate_files'].append(entry)
            else:
                record['retained_standard_perambulators'].append(entry)
    assert set(intermediates) == {BASE / f'cfg1000-t{t}-unsmeared.npy' for t in (2, 3)}
    opened = subprocess.run(['/usr/sbin/lsof', '-Fpn', '--', *map(str, intermediates)],
                            text=True, capture_output=True)
    assert opened.returncode == 1 and not opened.stdout, opened.stdout + opened.stderr
    record['intermediate_files_not_open'] = True
    for path in intermediates:
        path.unlink()
    record['status'] = 'reclaimed_verified_intermediates'
else:
    record['status'] = 'no_cleanup_needed'
record['free_disk_GiB_after'] = shutil.disk_usage(BASE).free / 2**30
(REPORT / 'intermediate-cleanup.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
