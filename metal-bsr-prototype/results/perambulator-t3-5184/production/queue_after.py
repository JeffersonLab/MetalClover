from pathlib import Path
import datetime
import hashlib
import json
import os
import shutil
import signal
import subprocess
import time

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
WORK = Path('/private/tmp/chroma-time3-production')
REPORT = ROOT / 'metal-bsr-prototype/results/perambulator-t3-5184'
BASE = Path('/Users/admin/Downloads/MetalClover-results')
PREVIOUS = BASE / 'cfg1000-t2-unsmeared'
POINT = BASE / 'cfg1000-t3-unsmeared'
STATUS = Path(str(POINT) + '.queue.json')
if STATUS.exists():
    raise RuntimeError('time-three queue already exists')


def read(path):
    return json.loads(Path(path).read_text())


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


record = {
    'status': 'queued',
    'source_time': 3,
    'sink_time': 3,
    'physical_rhs': 20736,
    'queue_PID': os.getpid(),
    'dependency_pipeline': str(PREVIOUS) + '.pipeline.json',
    'pipeline_script': str(WORK / 'full_run.py'),
    'queued_UTC': utc_now(),
    'continues_without_iPhone_connection': True,
    'minimum_free_disk_GiB_before_start': 17.0,
}


def save():
    temporary = Path(str(STATUS) + '.tmp')
    temporary.write_text(json.dumps(record, indent=2) + '\n')
    temporary.replace(STATUS)


stop_requested = []
child = None


def on_signal(sig, _):
    stop_requested.append(sig)
    if child is not None and child.poll() is None:
        child.send_signal(signal.SIGTERM)


signal.signal(signal.SIGTERM, on_signal)
signal.signal(signal.SIGINT, on_signal)
keeper = subprocess.Popen(
    ['/usr/bin/caffeinate', '-i', '-w', str(os.getpid())],
    stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
    stderr=subprocess.DEVNULL,
)
record['awake_keeper_PID'] = keeper.pid
save()
try:
    while not stop_requested:
        dependency = read(str(PREVIOUS) + '.pipeline.json')
        free_gib = shutil.disk_usage(BASE).free / 2**30
        record.update(
            dependency_status=dependency.get('status'),
            dependency_stage=dependency.get('stage'),
            free_disk_GiB=free_gib,
            last_checked_UTC=utc_now(),
        )
        if dependency.get('status') != 'complete':
            record['waiting_for'] = 'time-two pipeline completion'
        elif free_gib < record['minimum_free_disk_GiB_before_start']:
            record['waiting_for'] = 'disk space'
        else:
            record['waiting_for'] = 'input validation'
            save()
            break
        save()
        time.sleep(30)
    if stop_requested:
        record['status'] = 'paused'
        save()
        raise SystemExit(75)

    # Pipeline completion includes independent DataLib key validation.
    first = read(str(PREVIOUS) + '.complete.json')
    peram = Path(dependency['standard_perambulator'])
    projected = read(str(peram) + '.complete.json')
    assert first['status'] == projected['status'] == 'complete'
    assert first['source_time'] == first['sink_time'] == 2
    assert first['completed_rhs'] == 20736
    assert first['session_count'] == first['MG_hierarchy_builds_across_sessions'] == 1
    assert first['maximum_full_relative'] <= 1e-12
    assert first['maximum_qprop_relative'] <= 1e-12
    assert projected['source_time'] == projected['sink_time'] == 2
    assert projected['spin_records'] == 16
    assert projected['full_gamma5_relative_Frobenius'] < 1e-12
    assert peram.stat().st_size == projected['output_bytes']
    staged = read(REPORT / 'staged-run.json')
    assert sha(staged['eigenvector_file']) == staged['eigenvector_SHA256']
    assert sha(staged['gauge_file']) == staged['gauge_SHA256']
    for name, value in read(REPORT / 'frozen-files.json').items():
        if name.startswith('/'):
            assert sha(name) == value, name
    record.update(
        status='running', waiting_for=None,
        dependency_standard_perambulator_SHA256=projected['output_SHA256'],
        started_UTC=utc_now(),
    )
    with (WORK / 'pipeline-controller.log').open('x') as log:
        child = subprocess.Popen(
            ['python3', str(WORK / 'full_run.py')],
            stdin=subprocess.DEVNULL, stdout=log, stderr=subprocess.STDOUT,
        )
        record['pipeline_PID'] = child.pid
        save()
        code = child.wait()
    record.update(
        status='complete' if code == 0 else ('paused' if code == 75 else 'failed'),
        exit_code=code, finished_UTC=utc_now(),
    )
    save()
    raise SystemExit(code)
except Exception as error:
    record.update(status='failed', error=str(error))
    save()
    raise
finally:
    if keeper.poll() is None:
        keeper.terminate()
    keeper.wait()
