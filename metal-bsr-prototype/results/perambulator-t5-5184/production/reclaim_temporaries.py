from pathlib import Path
import datetime
import hashlib
import json
import shutil
import subprocess

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
REPORT = ROOT / 'metal-bsr-prototype/results/perambulator-t5-5184'
BASE = Path('/Users/admin/Downloads/MetalClover-results')
entries = [
    {'path': '/private/tmp/mg-null24.bridge', 'bytes': 2717909072,
     'sha256': 'be9479cc3b7597b8ff2bc717055cacecadfff8e03e6ca37bcb892e6f6f392be7',
     'archived_provenance': '/private/tmp/mg-null24-bridge-provenance.json',
     'kind': 'regenerable temporary gauge/clover benchmark export'},
    {'path': '/private/tmp/mg-lowmem12.bridge', 'bytes': 127402064,
     'sha256': 'a3a5acfb2f27f99697be88caff2c751bd1e9ff3a4420549baf22f7706e66ae05',
     'archived_provenance': 'metal-bsr-prototype/results/near-null-svd12/settings.json',
     'kind': 'regenerable temporary gauge/clover benchmark export'},
    {'path': '/private/tmp/chroma-distillation-spin/spins1.reference.bin', 'bytes': 127401984,
     'sha256': 'fda09f3fb6d33198a0fa2546a27f4a4d096e2da658270cb2f99ce3f25bfdf1f0',
     'archived_provenance': 'metal-bsr-prototype/results/distillation-spins12/artifacts.json',
     'kind': 'regenerable obsolete four-spin benchmark reference'},
    {'path': '/private/tmp/chroma-eloy-spin-batching/spins1.reference.bin', 'bytes': 127401984,
     'sha256': '7dfe4d530a8033d943c5c8dedfa4b612c0df8b1b01028a02bbb2282ab52bcdb5',
     'archived_provenance': 'metal-bsr-prototype/results/eloy-spin-batching12/artifacts.json',
     'kind': 'regenerable obsolete four-spin benchmark reference'},
]


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


assert all(Path(entry['path']).is_file() for entry in entries)
for entry in entries:
    path = Path(entry['path'])
    assert path.stat().st_size == entry['bytes'], path
    assert sha(path) == entry['sha256'], path
    entry['archived_hash_verified'] = True
opened = subprocess.run(['/usr/sbin/lsof', '-Fpn', '--', *[entry['path'] for entry in entries]],
                        text=True, capture_output=True)
assert opened.returncode == 1 and not opened.stdout, opened.stdout + opened.stderr
record = {
    'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'reason': 'space for requested queued time-five production, retaining all standard perambulators and original inputs',
    'free_disk_GiB_before': shutil.disk_usage(BASE).free / 2**30,
    'files_not_open': True,
    'files': entries,
    'reclaimed_bytes': sum(entry['bytes'] for entry in entries),
    'all_original_inputs_standard_perambulators_and_active_run_files_retained': True,
    'active_time_four_source_executable_library_and_scripts_unchanged': True,
}
for entry in entries:
    Path(entry['path']).unlink()
record['free_disk_GiB_after'] = shutil.disk_usage(BASE).free / 2**30
(REPORT / 'temporary-artifact-cleanup.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
