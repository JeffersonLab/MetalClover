from pathlib import Path
import datetime
import hashlib
import json
import shutil
import subprocess

REPORT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU/metal-bsr-prototype/results/perambulator-t3-5184')
BASE = Path('/Users/admin/Downloads/MetalClover-results')
proposal = json.loads((REPORT / 'disk-space-proposal.json').read_text())
expected = {str(BASE / f'cfg1000-t{t}-unsmeared.npy') for t in (0, 1)}
assert {entry['path'] for entry in proposal['proposed_intermediate_files']} == expected


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


for entry in proposal['proposed_intermediate_files'] + proposal['preserved_standard_perambulators']:
    path = Path(entry['path'])
    assert path.stat().st_size == entry['bytes'], path
    assert sha(path) == entry['accepted_SHA256'], path
    entry['current_SHA256_matches_accepted_result'] = True

opened = subprocess.run(['/usr/sbin/lsof', '-Fpn', '--', *sorted(expected)],
                        text=True, capture_output=True)
assert opened.returncode == 1 and not opened.stdout, opened.stdout + opened.stderr
record = {
    'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'authorization': 'User explicitly approved removal of these two intermediate .npy files, retaining every standard perambulator',
    'reason': 'space for queued time-three production',
    'free_disk_GiB_before': shutil.disk_usage(BASE).free / 2**30,
    'removed_intermediate_files': proposal['proposed_intermediate_files'],
    'retained_standard_perambulators': proposal['preserved_standard_perambulators'],
    'all_original_gauge_eigenvectors_logs_and_checksums_retained': True,
    'intermediate_files_not_open': True,
}
for entry in proposal['proposed_intermediate_files']:
    Path(entry['path']).unlink()
record['free_disk_GiB_after'] = shutil.disk_usage(BASE).free / 2**30
(REPORT / 'intermediate-cleanup.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
