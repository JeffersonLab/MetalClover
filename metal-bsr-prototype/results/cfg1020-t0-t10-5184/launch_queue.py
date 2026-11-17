"""Launch the ten already validated workers with a continuously checked disk guard."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import subprocess
import time

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
REPORT = ROOT / 'metal-bsr-prototype/results/cfg1020-t0-t10-5184'
WORK = Path('/private/tmp/chroma-cfg1020-production')
BASE = Path('/Users/admin/Downloads/MetalClover-results')


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def main():
    assert not (REPORT / 'queue-launch.json').exists()
    gate = read(WORK / 'batch-storage-ready.json')
    assert gate['status'] in ('pass', 'insufficient_space', 'waiting_for_stable_snapshot')
    assert gate['all_standard_perambulators_retained'] is True
    driver = WORK / 'batch_queue.py'
    assert driver.read_bytes() == (REPORT / driver.name).read_bytes()
    tested = read(REPORT / 'queue-safety-validation.json')
    assert tested['status'] == 'pass'
    assert hashlib.sha256(driver.read_bytes()).hexdigest() == tested['driver_SHA256']
    assert read(REPORT / 'frozen-runtime-recheck.json')['status'] == 'pass'
    for t in range(1, 11):
        assert not (BASE / f'cfg1020-t{t}-unsmeared.queue.json').exists()
        assert not Path(f'/private/tmp/chroma-cfg1020-time{t}-production/queue-controller.log').exists()
    record = {'status': 'launching', 'configuration': 1020,
              'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_sink_times': list(range(1, 11)),
              'physical_RHS_per_time': 20736,
              'additional_physical_RHS_total': 207360,
              'detached_local_workers': True,
              'storage_status_at_launch': gate['status'],
              'storage_pending_at_launch': not gate.get('forecast_all_runs_fit', False),
              'all_standard_perambulators_retained': True,
              'jobs': []}
    write(REPORT / 'queue-launch.json', record)
    children = []
    for t in range(1, 11):
        log = Path(f'/private/tmp/chroma-cfg1020-time{t}-production/queue-controller.log')
        with log.open('x') as stream:
            child = subprocess.Popen(['python3', str(driver), '--time', str(t)],
                                     cwd=ROOT, stdin=subprocess.DEVNULL,
                                     stdout=stream, stderr=subprocess.STDOUT,
                                     start_new_session=True)
        children.append(child)
        record['jobs'].append({'source_time': t, 'sink_time': t,
                               'queue_PID': child.pid, 'log': str(log)})
        write(REPORT / 'queue-launch.json', record)
    for attempt in range(30):
        assert all(child.poll() is None for child in children), 'Queue worker exited at launch'
        statuses = [BASE / f'cfg1020-t{t}-unsmeared.queue.json' for t in range(1, 11)]
        if all(path.exists() and read(path).get('batch_storage_status') in ('pass', 'insufficient_space', 'waiting_for_stable_snapshot')
               for path in statuses):
            break
        time.sleep(1)
    for job in record['jobs']:
        t = job['source_time']
        status = read(BASE / f'cfg1020-t{t}-unsmeared.queue.json')
        assert status['status'] == 'queued'
        assert status['queue_PID'] == job['queue_PID']
        assert status['source_time'] == status['sink_time'] == t
        assert status['batch_storage_status'] in ('pass', 'insufficient_space', 'waiting_for_stable_snapshot')
        os.kill(job['queue_PID'], 0)
        os.kill(status['awake_keeper_PID'], 0)
        job.update(awake_keeper_PID=status['awake_keeper_PID'],
                   dependency_configuration=1000 if t == 1 else 1020,
                   dependency_time=10 if t == 1 else t - 1, verified_status='queued')
        write(REPORT.parent / f'cfg1020-perambulator-t{t}-5184/queue-launch.json', job)
        write(REPORT.parent / f'cfg1020-perambulator-t{t}-5184/queue-initial.json', status)
    record['status'] = 'queued'
    write(REPORT / 'queue-launch.json', record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
