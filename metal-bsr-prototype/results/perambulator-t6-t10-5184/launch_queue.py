"""Launch the five already validated workers once the disk forecast passes."""
from pathlib import Path
import datetime
import hashlib
import json
import os
import subprocess
import time

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
REPORT = ROOT / 'metal-bsr-prototype/results/perambulator-t6-t10-5184'
WORK = Path('/private/tmp/chroma-times6-10-production')
BASE = Path('/Users/admin/Downloads/MetalClover-results')


def read(path):
    return json.loads(Path(path).read_text())


def write(path, value):
    Path(path).write_text(json.dumps(value, indent=2) + '\n')


def main():
    assert not (REPORT / 'queue-launch.json').exists()
    gate = read(WORK / 'batch-storage-ready.json')
    assert gate['status'] == 'pass' and gate['forecast_all_runs_fit']
    assert gate['minimum_new_job_start_free_GiB'] >= 18
    assert gate['minimum_forecast_free_GiB'] >= 4
    driver = WORK / 'batch_queue.py'
    assert driver.read_bytes() == (REPORT / driver.name).read_bytes()
    tested = read(REPORT / 'queue-safety-validation.json')
    assert tested['status'] == 'pass'
    assert hashlib.sha256(driver.read_bytes()).hexdigest() == tested['driver_SHA256']
    assert read(REPORT / 'frozen-runtime-recheck.json')['status'] == 'pass'
    for t in range(6, 11):
        assert not (BASE / f'cfg1000-t{t}-unsmeared.queue.json').exists()
        assert not Path(f'/private/tmp/chroma-time{t}-production/queue-controller.log').exists()
    record = {'status': 'launching',
              'UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'source_sink_times': list(range(6, 11)),
              'physical_RHS_per_time': 20736,
              'additional_physical_RHS_total': 103680,
              'detached_local_workers': True,
              'all_standard_perambulators_retained': True,
              'jobs': []}
    write(REPORT / 'queue-launch.json', record)
    children = []
    for t in range(6, 11):
        log = Path(f'/private/tmp/chroma-time{t}-production/queue-controller.log')
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
        statuses = [BASE / f'cfg1000-t{t}-unsmeared.queue.json' for t in range(6, 11)]
        if all(path.exists() and read(path).get('batch_storage_status') == 'pass'
               for path in statuses):
            break
        time.sleep(1)
    for job in record['jobs']:
        t = job['source_time']
        status = read(BASE / f'cfg1000-t{t}-unsmeared.queue.json')
        assert status['status'] == 'queued'
        assert status['queue_PID'] == job['queue_PID']
        assert status['source_time'] == status['sink_time'] == t
        assert status['batch_storage_status'] == 'pass'
        os.kill(job['queue_PID'], 0)
        os.kill(status['awake_keeper_PID'], 0)
        job.update(awake_keeper_PID=status['awake_keeper_PID'],
                   dependency_time=t - 1, verified_status='queued')
        write(REPORT.parent / f'perambulator-t{t}-5184/queue-launch.json', job)
        write(REPORT.parent / f'perambulator-t{t}-5184/queue-initial.json', status)
    record['status'] = 'queued'
    write(REPORT / 'queue-launch.json', record)
    print(json.dumps(record, indent=2))


if __name__ == '__main__':
    main()
