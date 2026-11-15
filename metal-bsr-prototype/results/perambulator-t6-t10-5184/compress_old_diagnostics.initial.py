"""Preserve old diagnostics byte-for-byte while reducing APFS disk allocation."""

from pathlib import Path
import datetime
import hashlib
import json
import os
import shutil
import stat
import subprocess
import time

REPORT = Path(__file__).resolve().parent
BASE = Path('/Users/admin/Downloads/MetalClover-results')
TOOL = Path('/private/tmp/metalclover-storage-tools/afsctool')
TARGET_GIB = 55.0
ROOTS = [
    Path('/Users/admin/Downloads/redstar_kanoni.devel/runs'),
    Path('/Users/admin/Downloads/macos-datalib-reorg/evidence'),
    Path('/Users/admin/Downloads/chromaform.devel-datalib-reorg/build'),
]
EXCLUDED_DIRS = {'.git', 'sources', 'src', 'docs', 'data', 'reference',
                 'test_inputs', 'run_inputs', 'chroma_inputs', 'weights'}
MACHO = {b'\xcf\xfa\xed\xfe', b'\xce\xfa\xed\xfe', b'\xca\xfe\xba\xbe'}
CUTOFF = time.time() - 24 * 3600
JOURNAL = REPORT / 'storage-preservation.jsonl'
STATUS = REPORT / 'storage-preservation.json'


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def snapshot(info):
    return (info.st_dev, info.st_ino, info.st_size, info.st_mtime_ns,
            info.st_ctime_ns, info.st_mode, info.st_uid, info.st_gid)


def journal(record):
    with JOURNAL.open('a') as stream:
        stream.write(json.dumps(record) + '\n')
        stream.flush()
        os.fsync(stream.fileno())


def opened_paths():
    result = subprocess.run(['/usr/sbin/lsof', '-Fn'], capture_output=True, text=True)
    if result.returncode not in (0, 1):
        raise RuntimeError('Cannot inspect open files: ' + result.stderr)
    return {line[1:] for line in result.stdout.splitlines() if line.startswith('n/')}


def classify(path, info):
    if path.suffix == '.log' and info.st_size >= 2**20:
        return 0
    if path.suffix == '.a' and info.st_size >= 2**20:
        with path.open('rb') as stream:
            return 1 if stream.read(8) == b'!<arch>\n' else None
    if not path.suffix and info.st_size >= 4 * 2**20:
        with path.open('rb') as stream:
            return 1 if stream.read(4) in MACHO else None
    if path.suffix == '.Po' and info.st_size >= 128 * 1024:
        return 2
    if path.suffix == '.o' and info.st_size >= 2**20:
        with path.open('rb') as stream:
            return 3 if stream.read(4) in MACHO else None
    return None


if JOURNAL.exists() or STATUS.exists():
    raise RuntimeError('Storage audit already exists; inspect it before another run')
assert TOOL.is_file()
candidates = []
for root in ROOTS:
    for parent, dirs, files in os.walk(root):
        dirs[:] = [name for name in dirs if name not in EXCLUDED_DIRS
                   and not name.startswith('source')]
        for name in files:
            path = Path(parent, name)
            if path.is_symlink() or '.metalclover-' in name:
                continue
            info = path.stat()
            if not stat.S_ISREG(info.st_mode) or info.st_nlink != 1:
                continue
            if info.st_flags or info.st_mtime > CUTOFF or info.st_blocks == 0:
                continue
            if info.st_size >= 2**31:
                continue
            category = classify(path, info)
            if category is not None:
                candidates.append((category, -info.st_size, path))
candidates.sort()
(REPORT / 'storage-candidates.json').write_text(json.dumps(
    [{'priority': category, 'path': str(path), 'bytes': -size}
     for category, size, path in candidates], indent=2) + '\n')

summary = {
    'status': 'running', 'started_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
    'method': 'APFS transparent compression on a temporary clone; complete SHA256 before/after; atomic replacement with original retained until final readback passes',
    'all_original_file_paths_and_byte_contents_preserved': True,
    'numerical_datasets_and_source_files_excluded': True,
    'files_compressed': 0, 'files_skipped': 0, 'logical_bytes_preserved': 0,
    'allocated_bytes_saved_per_file': 0, 'target_free_GiB': TARGET_GIB,
    'free_disk_GiB_before': shutil.disk_usage(BASE).free / 2**30,
    'candidate_count': len(candidates),
}


def save():
    summary['free_disk_GiB'] = shutil.disk_usage(BASE).free / 2**30
    temp = Path(str(STATUS) + '.tmp')
    temp.write_text(json.dumps(summary, indent=2) + '\n')
    temp.replace(STATUS)


opened = opened_paths()
checked = time.monotonic()
save()
for _, _, path in candidates:
    if shutil.disk_usage(BASE).free / 2**30 >= TARGET_GIB:
        break
    if time.monotonic() - checked > 30:
        opened = opened_paths()
        checked = time.monotonic()
    if str(path) in opened:
        summary['files_skipped'] += 1
        continue
    original = path.stat()
    if original.st_mtime > CUTOFF or original.st_nlink != 1 or original.st_flags:
        summary['files_skipped'] += 1
        continue
    temporary = path.with_name('.' + path.name + f'.metalclover-compressed-{os.getpid()}')
    backup = path.with_name('.' + path.name + f'.metalclover-original-{os.getpid()}')
    assert not temporary.exists() and not backup.exists()
    digest = sha(path)
    entry = {'path': str(path), 'bytes': original.st_size, 'SHA256': digest,
             'allocated_before': original.st_blocks * 512,
             'temporary': str(temporary), 'original_backup': str(backup)}
    journal(dict(entry, phase='prepared'))
    installed = False
    try:
        subprocess.run(['/bin/cp', '-c', '-p', str(path), str(temporary)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        result = subprocess.run(['/usr/bin/nice', '-n', '15', str(TOOL), '-c',
                                 '-T', 'ZLIB', '-1', '-s', '10', '-m', '2147483647',
                                 str(temporary)], capture_output=True, text=True)
        if result.returncode or not temporary.stat().st_flags & 32:
            journal(dict(entry, phase='skipped', reason='compression not applied',
                         tool_exit=result.returncode, tool_stderr=result.stderr[-1000:]))
            temporary.unlink()
            summary['files_skipped'] += 1
            continue
        os.utime(temporary, ns=(original.st_atime_ns, original.st_mtime_ns))
        after = temporary.stat()
        assert after.st_size == original.st_size and sha(temporary) == digest
        assert stat.S_IMODE(after.st_mode) == stat.S_IMODE(original.st_mode)
        assert (after.st_uid, after.st_gid) == (original.st_uid, original.st_gid)
        if snapshot(path.stat()) != snapshot(original):
            temporary.unlink()
            journal(dict(entry, phase='skipped', reason='original changed during preparation'))
            summary['files_skipped'] += 1
            continue
        if str(path) in opened_paths():
            temporary.unlink()
            journal(dict(entry, phase='skipped', reason='file was opened during preparation'))
            summary['files_skipped'] += 1
            continue
        journal(dict(entry, phase='verified_temporary', allocated_after=after.st_blocks * 512))
        os.link(path, backup)
        os.replace(temporary, path)
        installed = True
        assert sha(path) == digest and path.stat().st_size == original.st_size
        assert path.stat().st_flags & 32
        with path.open('rb') as stream:
            os.fsync(stream.fileno())
        backup.unlink()
        saved = original.st_blocks * 512 - path.stat().st_blocks * 512
        summary['files_compressed'] += 1
        summary['logical_bytes_preserved'] += original.st_size
        summary['allocated_bytes_saved_per_file'] += saved
        journal(dict(entry, phase='committed', allocated_after=path.stat().st_blocks * 512,
                     full_final_path_SHA256_matches=True,
                     mtime_and_permissions_preserved=True,
                     free_disk_GiB=shutil.disk_usage(BASE).free / 2**30))
        save()
        if summary['files_compressed'] % 10 == 0:
            print(json.dumps({'files': summary['files_compressed'],
                              'free_GiB': round(summary['free_disk_GiB'], 3)}), flush=True)
    except Exception as error:
        if backup.exists():
            if installed:
                os.replace(backup, path)
            else:
                backup.unlink()
        if temporary.exists():
            temporary.unlink()
        journal(dict(entry, phase='rolled_back', error=str(error)))
        summary.update(status='failed', error=str(error), failed_path=str(path))
        save()
        raise
    finally:
        if temporary.exists():
            temporary.unlink()

summary.update(status='complete', finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
save()
print(json.dumps(summary, indent=2), flush=True)
