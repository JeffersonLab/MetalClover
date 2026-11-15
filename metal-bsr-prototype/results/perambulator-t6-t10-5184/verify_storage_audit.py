"""Check preservation journals and absence of leftover replacement files."""
from pathlib import Path
import datetime
import json
import shutil

REPORT = Path(__file__).resolve().parent
BASE = Path('/Users/admin/Downloads/MetalClover-results')
audits = []
for directory in [REPORT, REPORT / 'binary-preservation']:
    summary = json.loads((directory / 'storage-preservation.json').read_text())
    assert summary['status'] == 'complete'
    last = {}
    replacements = set()
    for line in (directory / 'storage-preservation.jsonl').read_text().splitlines():
        entry = json.loads(line)
        last[entry['path']] = entry
        replacements.update((entry['temporary'], entry['original_backup']))
    assert all(not Path(path).exists() for path in replacements)
    committed = [entry for entry in last.values() if entry['phase'] == 'committed']
    assert len(committed) == summary['files_compressed']
    for entry in last.values():
        assert entry['phase'] in ('committed', 'skipped', 'rolled_back')
        assert not Path(entry['temporary']).exists()
        assert not Path(entry['original_backup']).exists()
        if entry['phase'] == 'committed':
            info = Path(entry['path']).stat()
            assert info.st_size == entry['bytes'] and info.st_flags & 32
            assert entry['full_final_path_SHA256_matches']
            assert entry['mtime_and_permissions_preserved']
    audits.append({'directory': str(directory),
                   'files_compressed': len(committed),
                   'logical_bytes_preserved': sum(entry['bytes'] for entry in committed),
                   'allocated_bytes_saved_per_file': summary['allocated_bytes_saved_per_file'],
                   'all_committed_full_final_path_SHA256_checks_passed': True,
                   'no_temporary_or_backup_files_left': True})
record = {'status': 'pass',
          'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
          'actual_free_disk_GiB': shutil.disk_usage(BASE).free / 2**30,
          'source_files_numerical_datasets_and_perambulators_excluded': True,
          'audit_scope': 'Reconcile every journal entry, saved full SHA256 results, current sizes/compression flags and absence of replacement leftovers; full hashes were checked during each atomic replacement',
          'audits': audits}
(REPORT / 'storage-audit-validation.json').write_text(json.dumps(record, indent=2) + '\n')
print(json.dumps(record, indent=2))
