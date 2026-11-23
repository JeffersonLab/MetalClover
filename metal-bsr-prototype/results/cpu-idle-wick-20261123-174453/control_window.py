"""Suspend only verified queue processes; always resume them within one hour."""
from pathlib import Path
from datetime import datetime, timezone
import hashlib
import json
import os
import signal
import subprocess
import time

WORK = Path(__file__).resolve().parent
PLAN = json.loads((WORK / 'plan.json').read_text())
REPORT = Path(PLAN['report_directory'])


def utc():
    return datetime.now(timezone.utc).isoformat()


def process(pid):
    result = subprocess.run(['/bin/ps', '-p', str(pid), '-o',
                             'pid=,ppid=,stat=,lstart=,command='],
                            capture_output=True, text=True)
    if result.returncode:
        return None
    fields = result.stdout.strip().split(None, 8)
    if len(fields) != 9:
        raise RuntimeError('Cannot identify process ' + str(pid))
    return {'pid': int(fields[0]), 'parent': int(fields[1]), 'state': fields[2],
            'started': ' '.join(fields[3:8]), 'command': fields[8],
            'command_SHA256': hashlib.sha256(fields[8].encode()).hexdigest()}


def same_process(expected, actual):
    return actual is not None and all(expected[key] == actual[key]
                                      for key in ['pid', 'started', 'command_SHA256'])


def read(path):
    try:
        return json.loads(path.read_text())
    except FileNotFoundError:
        return {}


record = {'window_id': PLAN['window_id'], 'status': 'validating',
          'controller_PID': os.getpid(), 'UTC': utc(), 'targets': [],
          'method': PLAN['method'], 'one_hierarchy_and_live_solve_preserved': True}


def save():
    record['updated_UTC'] = utc()
    for path in [WORK / 'status.json', REPORT / 'pause-status.json', Path(PLAN['global_status_file'])]:
        temporary = path.with_suffix('.tmp')
        with temporary.open('w') as stream:
            stream.write(json.dumps(record, indent=2) + '\n')
            stream.flush()
            os.fsync(stream.fileno())
        temporary.replace(path)


stopped = []
requested = []
for sig in [signal.SIGINT, signal.SIGTERM]:
    signal.signal(sig, lambda number, frame: requested.append(number))
keeper = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(os.getpid())],
                          stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                          stderr=subprocess.DEVNULL)
reason = 'setup_failure'
try:
    # Validate the entire set before sending any stop signal.
    for target in PLAN['targets']:
        current = process(target['pid'])
        assert current is not None and not current['state'].startswith(('T', 'Z')), target
        assert all(part in current['command'] for part in target['command_contains']), target
        item = dict(target, original_identity=current, stopped=False, resumed=False)
        record['targets'].append(item)
    save()
    for item in record['targets']:
        assert same_process(item['original_identity'], process(item['pid']))
        os.kill(item['pid'], signal.SIGSTOP)
        stopped.append(item)
        for attempt in range(20):
            current = process(item['pid'])
            if current is not None and current['state'].startswith('T'):
                break
            time.sleep(0.05)
        assert current is not None and current['state'].startswith('T'), item
        item['stopped'] = True
    started = time.monotonic()
    deadline = started + PLAN['maximum_pause_seconds']
    record.update(status='paused', paused_UTC=utc(),
                  maximum_pause_seconds=PLAN['maximum_pause_seconds'],
                  automatic_resume_deadline_UTC=datetime.fromtimestamp(
                      time.time() + PLAN['maximum_pause_seconds'], timezone.utc).isoformat())
    save()
    print(json.dumps({'status': 'paused', 'window_id': PLAN['window_id'],
                      'processes': len(stopped),
                      'automatic_resume_deadline_UTC': record['automatic_resume_deadline_UTC']}), flush=True)
    reason = 'one_hour_failsafe'
    while time.monotonic() < deadline:
        if requested:
            reason = 'pause_controller_stop_requested'
            break
        result = read(WORK / 'resume-request.json')
        if result.get('window_id') == PLAN['window_id']:
            reason = result.get('reason', 'explicit_resume_request')
            break
        time.sleep(5)
    record['pause_seconds'] = time.monotonic() - started
except BaseException as error:
    record['error'] = str(error)
    reason = 'pause_controller_error'
    raise
finally:
    record.update(status='resuming', resume_reason=reason, resume_started_UTC=utc())
    save()
    errors = []
    # Reverse pause order: native solve, resource monitor, parents, then gates.
    for item in reversed(stopped):
        try:
            current = process(item['pid'])
            assert same_process(item['original_identity'], current), 'Process identity changed'
            os.kill(item['pid'], signal.SIGCONT)
            item['resumed'] = True
        except Exception as error:
            errors.append({'role': item['role'], 'pid': item['pid'], 'error': str(error)})
    record.update(status='resumed' if not errors else 'needs_attention',
                  resumed_UTC=utc(), resume_errors=errors)
    save()
    keeper.terminate()
    keeper.wait()
    print(json.dumps({'status': record['status'], 'resume_reason': reason,
                      'pause_seconds': record.get('pause_seconds'), 'resume_errors': errors}), flush=True)
