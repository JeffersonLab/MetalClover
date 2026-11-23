"""Recover completed time6 and resume the original, unchanged numerical jobs."""
from pathlib import Path
import datetime
import hashlib
import importlib.util
import json
import os
import shutil
import subprocess
import sys
import time

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
RESULTS = ROOT / 'metal-bsr-prototype/results'
HERE = Path(__file__).resolve().parent
BATCH = RESULTS / 'cfg1020-t0-t10-5184'
WORK = Path('/private/tmp/chroma-cfg1020-production')
BASE = Path('/Users/admin/Downloads/MetalClover-results')
BACKUP = Path('/Users/admin/Downloads/MetalClover-runtime/cfg1020-frozen-20261123')
POINT = BASE / 'cfg1020-t6-unsmeared'


def now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    path = Path(path)
    tmp = path.with_name(path.name + f'.tmp-{os.getpid()}')
    with tmp.open('w') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')
        stream.flush()
        os.fsync(stream.fileno())
    tmp.replace(path)


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def restore(source, target, expected, mode=None):
    source, target = Path(source), Path(target)
    assert sha(source) == expected, source
    if target.exists():
        assert sha(target) == expected, target
        return False
    target.parent.mkdir(parents=True, exist_ok=True)
    tmp = target.with_name(target.name + f'.restore-{os.getpid()}')
    shutil.copyfile(source, tmp)
    tmp.chmod(mode if mode is not None else source.stat().st_mode & 0o777)
    assert sha(tmp) == expected
    # Exclusive hard-link publication never replaces an unexpected live file.
    os.link(tmp, target)
    tmp.unlink()
    return True


def prepare():
    assert not BACKUP.exists(), 'Durable backup already prepared; inspect before reusing'
    sources = {x['target']: (Path(x['source']), x['mode'])
               for x in read(HERE / 'exact-restoration-plan.json')['restorations']}
    expected = {}
    for t in range(1, 11):
        report = RESULTS / f'cfg1020-perambulator-t{t}-5184'
        local = Path(f'/private/tmp/chroma-cfg1020-time{t}-production')
        expected.update(read(report / 'frozen-files.json'))
        for name in ['action.xml', 'full_run.py', 'capture_completed.py', 'render_completed.py',
                     'frozen-files.json', 'staged-run.json']:
            source = report / (('production/' + name) if name.endswith('.py') else name)
            target = str(local / name)
            sources[target] = (source, 0o644)
            if target not in expected:
                expected[target] = sha(source)
    for name in ['batch_queue.py', 'storage_watch.py', 'batch_status.py']:
        target = str(WORK / name)
        sources[target] = (BATCH / name, 0o644)
        if target not in expected:
            expected[target] = sha(BATCH / name)
    inventory = []
    for target, digest in expected.items():
        if not target.startswith('/'):
            continue
        p = Path(target)
        existed = p.exists()
        if existed:
            assert sha(p) == digest, target
        else:
            source, mode = sources[target]
            restore(source, p, digest, mode)
        backup = BACKUP / target.lstrip('/')
        backup.parent.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(p, backup)
        backup.chmod(p.stat().st_mode & 0o777)
        assert sha(backup) == digest
        inventory.append({'path': target, 'SHA256': digest, 'backup': str(backup),
                          'mode': p.stat().st_mode & 0o777, 'restored': not existed})
        os.utime(p, None)
    write(BACKUP / 'manifest.json', inventory)
    write(HERE / 'restoration-20261123.json', {'status': 'pass', 'UTC': now(),
          'restored_files': sum(x['restored'] for x in inventory),
          'all_frozen_bytes_unchanged': True, 'durable_backup': str(BACKUP), 'files': inventory})
    print('Restored and backed up', len(inventory), 'exact files', flush=True)


def supervise():
    manifest = read(BACKUP / 'manifest.json')
    storage = module('original_storage', WORK / 'storage_watch.py')
    while True:
        repaired = []
        for item in manifest:
            p = Path(item['path'])
            if restore(item['backup'], p, item['SHA256'], item['mode']):
                repaired.append(str(p))
            os.utime(p, None)
            for parent in [p.parent, p.parent.parent]:
                os.utime(parent, None)
        sample = storage.sample()
        sample['supervisor_PID'] = os.getpid()
        # Full-queue forecast remains honest. Each worker additionally receives
        # a scoped forecast through its own time, so a later disk shortage need
        # not block work whose complete output and reserve already fit.
        for name in ['batch-storage-ready.json', 'disk-budget-current.json', 'storage-supervisor.json']:
            write(WORK / name, sample)
        for t in range(7, 11):
            if 'snapshot' not in sample:
                plan = dict(sample)
            else:
                scope = [s for s in sample['snapshot']
                         if s['configuration'] != 1020 or s['time'] <= t]
                plan = storage.forecast(sample['actual_free_disk_GiB'] * 2**30, scope)
                plan.update(UTC=now(), supervisor_PID=os.getpid(),
                            forecast_scope=f'All remaining diagonal jobs through cfg1020 time{t}',
                            full_queue_forecast_all_runs_fit=sample['forecast_all_runs_fit'])
            write(WORK / f'storage-through-time{t}.json', plan)
        write(HERE / 'runtime-supervisor.json', {'status': 'running', 'PID': os.getpid(),
              'UTC': now(), 'repaired_missing_files_this_check': repaired,
              'all_backup_and_runtime_hashes_verified': True,
              'full_queue_additional_free_space_GiB': sample.get('additional_free_space_required_GiB')})
        if all(read(BASE / f'cfg1020-t{t}-unsmeared.queue.json').get('status') == 'complete'
               for t in range(7, 11)):
            write(HERE / 'runtime-supervisor.json', {'status': 'complete', 'UTC': now()})
            return
        time.sleep(30)


def launch():
    record_path = HERE / 'resume-launch-20261123.json'
    assert not record_path.exists(), 'Resume has already been launched'
    record = {'status': 'launching', 'UTC': now(), 'user_authorized_resume': True,
              'continue_without_chat_or_phone': True, 'configuration': 1020,
              'scope': 'Recover diagonal6; resume diagonal7..10; no source-zero or cfg1000 work',
              'durable_runtime_backup': str(BACKUP), 'jobs': []}
    # Archive pause records, and archive (rather than discard) the stopped workers.
    for marker in [BASE / 'cfg1020.USER_PAUSED.json', WORK / 'user-pause.json']:
        if marker.exists():
            shutil.copyfile(marker, HERE / ('archived-20261123-' + marker.name))
            value = read(marker)
            value.update(status='resumed', resumed_UTC=now(),
                         automatic_resume_authorized=True, resumed_on_explicit_user_request=True)
            write(marker, value)
    for t in range(7, 11):
        prefix = BASE / f'cfg1020-t{t}-unsmeared'
        assert not Path(str(prefix) + '.npy').exists()
        assert not Path(str(prefix) + '.pipeline.json').exists()
        status = Path(str(prefix) + '.queue.json')
        assert read(status)['status'] == 'paused'
        status.rename(status.with_name(status.name + '.paused-20261121'))
    write(record_path, record)
    for name, args in [('runtime-supervisor', ['supervise']), ('time6-postprocessing', ['recover6'])] + [
            (f'time{t}-queue', ['queue', str(t)]) for t in range(7, 11)]:
        log = HERE / f'{name}-20261123.log'
        with log.open('x') as stream:
            child = subprocess.Popen(['/usr/bin/caffeinate', '-i', 'python3', str(Path(__file__).resolve())] + args,
                                     stdin=subprocess.DEVNULL, stdout=stream, stderr=subprocess.STDOUT,
                                     cwd=ROOT, start_new_session=True)
        record['jobs'].append({'name': name, 'PID': child.pid, 'log': str(log)})
        write(record_path, record)
    record.update(status='launched', UTC=now())
    write(record_path, record)
    (WORK / 'ACTIVE.md').write_text(
        'RESUMED on explicit user request 2026-11-23. Detached local jobs continue without phone/chat.\n'
        'Scope: finish cfg1020 diagonal6 postprocessing WITHOUT reinversion; diagonal7..10 follow.\n'
        'All original FP64 checks and standard perambulators retained. No cfg1000 11..20 or P(t,0).\n'
        'Durable exact runtime backup: ' + str(BACKUP) + '\n'
        'Runtime guardian repairs only missing hash-verified files and refreshes temporary-file dates.\n'
        'Each worker uses storage-through-timeN.json; full-sequence forecast remains in storage-supervisor.json.\n'
        'If later jobs do not fit, they wait for space while the safe prefix proceeds.\n'
        'Recovery and live process records: ' + str(HERE) + '\n')
    print(json.dumps(record, indent=2))


def recover6():
    start = time.monotonic()
    archive = HERE / 'time6-before-recovery-20261123'
    archive.mkdir(exist_ok=False)
    for suffix in ['.queue.json', '.pipeline.json', '.run.json', '.session-000.status.json',
                   '.session-000.verify.log', '.state']:
        shutil.copyfile(str(POINT) + suffix, archive / ('point' + suffix))
    local = Path('/private/tmp/chroma-cfg1020-time6-production')
    for name in ['pipeline-controller.log', 'full-solve-controller.log']:
        shutil.copyfile(local / name, archive / name)
    run = read(str(POINT) + '.run.json')
    assert run['completed_rhs'] == 20736 and run['exit_code'] == 0
    assert len(list(BASE.glob(POINT.name + '.session-*.status.json'))) == 1
    assert not Path(str(POINT) + '.complete.json').exists()
    q = read(str(POINT) + '.queue.json')
    q.update(status='recovering', stage='postprocessing', recovery_PID=os.getpid(),
             original_queue_PID=q['queue_PID'], queue_PID=os.getpid(),
             original_failure_preserved=str(archive))
    for key in ['error', 'exit_code', 'finished_UTC', 'controller_PID', 'awake_keeper_PID']:
        q.pop(key, None)
    write(str(POINT) + '.queue.json', q)
    pipeline = read(str(POINT) + '.pipeline.json')
    pipeline.update(status='running', stage='recovery-verify-only', recovery_started_UTC=now(),
                    original_failure_preserved=str(archive),
                    recovery_does_not_repeat_inversions=True)
    write(str(POINT) + '.pipeline.json', pipeline)
    env = dict(os.environ, OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    for path, digest in run['inputs_sha256'].items():
        assert sha(path) == digest, path
    command = run['command'] + ['--resume']
    verify = Path(str(POINT) + '.session-000.verify.log')
    assert verify.stat().st_size == 0
    with verify.open('w') as stream:
        subprocess.run(command, stdin=subprocess.DEVNULL, stdout=stream,
                       stderr=subprocess.STDOUT, env=env, check=True)
    verified = verify.read_text()
    assert 'RESUME_VERIFIED completed_rhs=20736' in verified
    assert 'NEAR_NULL_SVD ' not in verified and 'MetalMGClover array: rhs=' not in verified
    controller = module('frozen_unsmeared',
                        Path('/private/tmp/chroma-times6-10-production/runtime/controllers/run_unsmeared.py'))
    result = controller.summarize(POINT)
    assert result['session_count'] == result['MG_hierarchy_builds_across_sessions'] == 1
    run.update(status='complete', payload_sha256=result['payload_sha256'],
               postverification_recovered_UTC=now(),
               recovery_preserved_original_solve_session=True)
    write(str(POINT) + '.run.json', run)
    for stage in ['project', 'reader']:
        pipeline.update(status='running', stage=stage)
        write(str(POINT) + '.pipeline.json', pipeline)
        with (local / f'full-{stage}-controller.log').open('x') as stream:
            subprocess.run(pipeline['commands'][stage], stdin=subprocess.DEVNULL,
                           stdout=stream, stderr=subprocess.STDOUT, check=True)
    pipeline.update(status='complete', stage='complete', finished_UTC=now(), exit_code=0,
                    recovery_active_seconds=time.monotonic() - start,
                    elapsed_seconds=(datetime.datetime.now(datetime.timezone.utc) -
                                     datetime.datetime.fromisoformat(pipeline['started_UTC'])).total_seconds(),
                    elapsed_seconds_includes_runtime_failure_and_user_pause=True,
                    numerical_solve_seconds=run['elapsed_seconds'])
    write(str(POINT) + '.pipeline.json', pipeline)
    subprocess.run(['python3', str(local / 'capture_completed.py'), '6'], check=True)
    batch = module('original_batch', WORK / 'batch_queue.py')
    batch.retire_intermediates(6, [6], 'intermediate-retirement.json')
    subprocess.run(['python3', str(local / 'render_completed.py'), '6'], check=True)
    q = read(str(POINT) + '.queue.json')
    q.update(status='complete', stage='complete', exit_code=0, finished_UTC=now(),
             recovered_after_original_pipeline_failure=True, recovery_PID=os.getpid(),
             recovery_active_seconds=time.monotonic() - start,
             complete_evidence_saved=True, point_intermediate_retired=True)
    write(RESULTS / 'cfg1020-perambulator-t6-5184/queue-complete.json', q)
    write(str(POINT) + '.queue.json', q)
    write(HERE / 'time6-recovered.json', {'status': 'pass', 'UTC': now(),
          'one_original_solve_session_and_hierarchy': True, 'no_repeated_inversions': True,
          'recovery_active_seconds': time.monotonic() - start,
          'maximum_full_relative': result['maximum_full_relative']})


def queue(t):
    batch = module('original_batch', WORK / 'batch_queue.py')
    batch.STORAGE_PLAN = WORK / f'storage-through-time{t}.json'
    batch.main(t)


if __name__ == '__main__':
    try:
        {'prepare': prepare, 'supervise': supervise, 'recover6': recover6, 'launch': launch,
         'queue': lambda: queue(int(sys.argv[2]))}[sys.argv[1]]()
    except Exception as error:
        write(HERE / (sys.argv[1] + '-failure-20261123.json'),
              {'status': 'failed', 'UTC': now(), 'error': repr(error), 'PID': os.getpid()})
        if sys.argv[1] == 'recover6':
            path = Path(str(POINT) + '.pipeline.json')
            value = read(path)
            value.update(status='failed', recovery_error=repr(error), recovery_failed_UTC=now())
            write(path, value)
        raise
