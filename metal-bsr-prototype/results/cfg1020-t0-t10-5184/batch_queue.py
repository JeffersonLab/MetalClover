"""Detached sequential perambulator jobs, with full acceptance before handoff."""

from pathlib import Path
import argparse
import datetime
import hashlib
import json
import math
import os
import shutil
import signal
import subprocess
import time

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
RESULTS = ROOT / 'metal-bsr-prototype/results'
BASE = Path('/Users/admin/Downloads/MetalClover-results')
POINT_BYTES = 6879707264
STANDARD_BYTES = 6907158528
NUMERICAL_PAYLOAD_BYTES = 16 * 5184**2 * 16
RHS = 20736
TOLERANCE = 1e-12
MINIMUM_FREE_GIB = 18.0
STORAGE_PLAN = Path('/private/tmp/chroma-cfg1020-production/batch-storage-ready.json')


def utc_now():
    return datetime.datetime.now(datetime.timezone.utc).isoformat()


def read_json(path):
    try:
        return json.loads(Path(path).read_text())
    except FileNotFoundError:
        return {}


def write_json(path, value):
    path = Path(path)
    temporary = path.with_name(path.name + f'.tmp-{os.getpid()}')
    with temporary.open('w') as stream:
        stream.write(json.dumps(value, indent=2) + '\n')
        stream.flush()
        os.fsync(stream.fileno())
    temporary.replace(path)


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


def paths(t):
    assert 1 <= t <= 10
    point = BASE / f'cfg1020-t{t}-unsmeared'
    standard = BASE / f'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb1020'
    report = RESULTS / f'cfg1020-perambulator-t{t}-5184'
    return point, standard, report


def dependency_ready(pipeline, queue):
    return pipeline.get('status') == queue.get('status') == 'complete'


def storage_ready(plan):
    return (plan.get('status') == 'pass'
            and plan.get('forecast_all_runs_fit') is True
            and plan.get('all_standard_perambulators_retained') is True
            and plan.get('minimum_forecast_free_GiB', 0) >= 4.0
            and plan.get('minimum_new_job_start_free_GiB', 0) >= MINIMUM_FREE_GIB)


def bounded(value):
    return math.isfinite(value) and 0 <= value <= TOLERANCE


def validate_completion(t):
    point, standard, _ = paths(t)
    accepted = read_json(str(point) + '.complete.json')
    projected = read_json(str(standard) + '.complete.json')
    pipeline = read_json(str(point) + '.pipeline.json')
    assert accepted.get('status') == projected.get('status') == pipeline.get('status') == 'complete'
    assert accepted['source_time'] == accepted['sink_time'] == t
    assert projected['source_time'] == projected['sink_time'] == t
    assert pipeline['source_time'] == pipeline['sink_time'] == t
    assert accepted['physical_rhs'] == accepted['completed_rhs'] == RHS
    assert accepted['num_vectors'] == 5184 and accepted['source_spins'] == 4
    assert accepted['session_count'] == accepted['MG_hierarchy_builds_across_sessions'] == 1
    assert bounded(accepted['maximum_full_relative'])
    assert bounded(accepted['maximum_qprop_relative'])
    assert accepted['full_CPU_FP64_tolerance'] == TOLERANCE
    assert accepted['payload_bytes'] == POINT_BYTES
    assert projected['spin_records'] == 16 and projected['num_vectors'] == 5184
    assert bounded(projected['full_gamma5_relative_Frobenius'])
    assert bounded(projected['max_QDP_dot_scaled_error'])
    # Native readback checks all matrices; metadata/page padding may vary by time.
    assert NUMERICAL_PAYLOAD_BYTES <= projected['output_bytes'] <= NUMERICAL_PAYLOAD_BYTES + 64 * 2**20
    assert standard.is_file() and not standard.is_symlink()
    assert standard.stat().st_size == projected['output_bytes']
    return accepted, projected


def verify_proof(t, accepted, projected):
    _, _, report = paths(t)
    proof = read_json(report / 'production-validation.json')
    assert proof.get('status') == 'pass' and proof.get('configuration') == 1020
    assert proof['source_time'] == proof['sink_time'] == t
    assert proof['physical_rhs'] == RHS and proof['source_vectors'] == 5184
    assert proof['hierarchies'] == proof['sessions'] == 1
    assert proof['native_array_calls'] == 1296 and proof['RHS_per_call'] == 16
    assert proof['point_SHA256'] == accepted['payload_sha256']
    assert proof['standard_perambulator_SHA256'] == projected['output_SHA256']
    assert proof['all_point_column_checksums_verified']
    assert proof['all_standard_record_values_readback_bitwise']
    assert proof['independent_DataLib_full_file_key_check'] == 'pass'
    assert proof['unchanged_frozen_source_executable_library_hashes']
    return proof


def retire_intermediates(owner, times, record_name):
    """Remove only accepted point intermediates; every standard file stays."""
    assert 1 <= owner <= 10 and all(1 <= t <= owner for t in times)
    _, _, destination = paths(owner)
    manifest = {'status': 'validating', 'UTC': utc_now(), 'owner_time': owner, 'configuration': 1020,
                'authorization': 'Requested unattended standard-perambulator sequence; retire only fully validated point intermediates after preserving their complete evidence',
                'all_standard_perambulators_and_original_inputs_retained': True,
                'free_disk_GiB_before': shutil.disk_usage(BASE).free / 2**30,
                'files': [], 'already_absent_times': []}
    planned = []
    for t in times:
        point, standard, report = paths(t)
        payload = Path(str(point) + '.npy')
        if not payload.exists():
            manifest['already_absent_times'].append(t)
            continue
        assert not payload.is_symlink()
        accepted, projected = validate_completion(t)
        proof = verify_proof(t, accepted, projected)
        before = payload.stat()
        assert before.st_size == POINT_BYTES
        assert sha(payload) == accepted['payload_sha256']
        assert sha(standard) == projected['output_SHA256']
        after = payload.stat()
        assert (before.st_ino, before.st_size, before.st_mtime_ns) == (
            after.st_ino, after.st_size, after.st_mtime_ns)
        opened = subprocess.run(['/usr/sbin/lsof', '-nP', '-Fpn', '--', str(payload)],
                                text=True, capture_output=True)
        assert opened.returncode == 1 and not opened.stdout, opened.stdout + opened.stderr
        entry = {'time': t, 'point_file': str(payload), 'point_bytes': POINT_BYTES,
                 'point_SHA256': accepted['payload_sha256'],
                 'retained_standard': str(standard), 'standard_bytes': projected['output_bytes'],
                 'standard_SHA256': projected['output_SHA256'],
                 'current_point_and_standard_hashes_verified': True,
                 'point_file_not_open': True, 'retired': False}
        manifest['files'].append(entry)
        planned.append((payload, report, proof, entry))
    # Validate the entire plan before removing any file.
    write_json(destination / record_name, manifest)
    for payload, report, proof, entry in planned:
        payload.unlink()
        entry['retired'] = True
        proof.update(point_intermediate_retired_after_full_validation=True,
                     retirement_record=f'../cfg1020-perambulator-t{owner}-5184/{record_name}',
                     standard_perambulator_retained_SHA256_rechecked=True)
        write_json(report / 'production-validation.json', proof)
        write_json(destination / record_name, manifest)
    manifest.update(status='complete', free_disk_GiB_after=shutil.disk_usage(BASE).free / 2**30)
    write_json(destination / record_name, manifest)
    return [entry['time'] for entry in manifest['files']]


class StageError(RuntimeError):
    def __init__(self, stage, code):
        super().__init__(f'{stage} exited with code {code}')
        self.code = code


def main(t):
    assert 1 <= t <= 10
    work = Path(f'/private/tmp/chroma-cfg1020-time{t}-production')
    point, _, report = paths(t)
    previous = BASE / 'cfg1000-t10-unsmeared' if t == 1 else paths(t - 1)[0]
    status_path = Path(str(point) + '.queue.json')
    if status_path.exists():
        raise RuntimeError(f'time-{t} queue already exists')
    record = {'status': 'queued', 'source_time': t, 'sink_time': t,
              'configuration': 1020, 'physical_rhs': RHS, 'queue_PID': os.getpid(), 'queued_UTC': utc_now(),
              'dependency_pipeline': str(previous) + '.pipeline.json',
              'pipeline_script': str(work / 'full_run.py'),
              'continues_without_iPhone_connection': True,
              'minimum_free_disk_GiB_before_start': MINIMUM_FREE_GIB,
              'finalize_own_evidence_before_success': True}
    child = None
    requested = []

    def save():
        write_json(status_path, record)

    def on_signal(sig, _):
        requested.append(sig)
        if child is not None and child.poll() is None:
            child.send_signal(signal.SIGTERM)

    def check_stop():
        if requested:
            raise StageError('explicit stop request', 75)

    def execute(stage, command):
        nonlocal child
        check_stop()
        record.update(stage=stage)
        save()
        with (work / f'{stage}-controller.log').open('x') as log:
            child = subprocess.Popen([str(x) for x in command], stdin=subprocess.DEVNULL,
                                     stdout=log, stderr=subprocess.STDOUT)
            record['controller_PID'] = child.pid
            if stage == 'pipeline':
                record['pipeline_PID'] = child.pid
            save()
            code = child.wait()
        if code:
            raise StageError(stage, code)

    def ensure_evidence(source_time):
        accepted, projected = validate_completion(source_time)
        _, _, source_report = paths(source_time)
        if read_json(source_report / 'production-validation.json').get('status') != 'pass':
            execute(f'capture-time{source_time}',
                    ['python3', work / 'capture_completed.py', source_time])
        verify_proof(source_time, accepted, projected)
        return accepted, projected

    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)
    keeper = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(os.getpid())],
                              stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                              stderr=subprocess.DEVNULL)
    record['awake_keeper_PID'] = keeper.pid
    save()
    started = None
    try:
        while True:
            check_stop()
            dependency = read_json(str(previous) + '.pipeline.json')
            dependency_queue = read_json(str(previous) + '.queue.json')
            storage = read_json(STORAGE_PLAN)
            record.update(dependency_status=dependency.get('status', 'not_started'),
                          dependency_stage=dependency.get('stage', 'not_started'),
                          dependency_queue_status=dependency_queue.get('status', 'not_started'),
                          batch_storage_status=storage.get('status', 'not_ready'),
                          waiting_for=('cfg1000 time10' if t == 1 else f'cfg1020 time{t-1}') + ' complete validation/finalization and disk budget',
                          last_checked_UTC=utc_now(), free_disk_GiB=shutil.disk_usage(BASE).free / 2**30)
            save()
            if dependency_ready(dependency, dependency_queue) and storage_ready(storage):
                break
            time.sleep(30)
        staged = read_json(work / 'staged-run.json')
        assert staged['source_time'] == staged['sink_time'] == t
        assert sha(staged['eigenvector_file']) == staged['eigenvector_SHA256']
        assert sha(staged['gauge_file']) == staged['gauge_SHA256']
        for name, expected in read_json(work / 'frozen-files.json').items():
            if name.startswith('/'):
                assert sha(name) == expected, name
        if t == 1:
            old_report = RESULTS / 'perambulator-t10-5184'
            accepted = read_json(str(previous) + '.complete.json')
            proof = read_json(old_report / 'production-validation.json')
            old_standard = BASE / 'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_10.sdb1000'
            projected = read_json(str(old_standard) + '.complete.json')
            assert accepted['status'] == projected['status'] == 'complete'
            assert accepted['source_time'] == accepted['sink_time'] == 10
            assert accepted['completed_rhs'] == accepted['physical_rhs'] == RHS
            assert accepted['session_count'] == accepted['MG_hierarchy_builds_across_sessions'] == 1
            assert bounded(accepted['maximum_full_relative']) and bounded(accepted['maximum_qprop_relative'])
            assert proof['status'] == 'pass' and proof['all_standard_record_values_readback_bitwise']
            assert proof['independent_DataLib_full_file_key_check'] == 'pass'
            assert projected['spin_records'] == 16 and bounded(projected['full_gamma5_relative_Frobenius'])
            assert old_standard.stat().st_size == projected['output_bytes']
            assert sha(old_standard) == projected['output_SHA256'] == proof['standard_perambulator_SHA256']
            write_json(report / 'cfg1000-time10-dependency.json',
                       {'status': 'pass', 'standard_SHA256': projected['output_SHA256'],
                        'all_original_full_FP64_checks_and_independent_reader_passed': True})
        else:
            _, projected = ensure_evidence(t - 1)
            shutil.copyfile(Path(str(previous) + '.queue.json'),
                            paths(t - 1)[2] / 'queue-complete.json')
        retired = retire_intermediates(t, range(1, t), 'dependency-retirement.json')
        for source_time in sorted(set(retired + ([t - 1] if t > 1 else []))):
            execute(f'render-dependency-time{source_time}',
                    ['python3', work / 'render_completed.py', source_time])
        while shutil.disk_usage(BASE).free / 2**30 < MINIMUM_FREE_GIB:
            check_stop()
            record.update(waiting_for='disk space', last_checked_UTC=utc_now(),
                          free_disk_GiB=shutil.disk_usage(BASE).free / 2**30)
            save()
            time.sleep(30)
        check_stop()
        started = time.monotonic()
        record.update(status='running', waiting_for=None, started_UTC=utc_now(),
                      dependency_standard_perambulator_SHA256=projected['output_SHA256'])
        save()
        execute('pipeline', ['python3', work / 'full_run.py'])
        accepted, projected = ensure_evidence(t)
        retire_intermediates(t, [t], 'intermediate-retirement.json')
        execute('render-completed', ['python3', work / 'render_completed.py', t])
        record.update(status='complete', stage='complete', finished_UTC=utc_now(),
                      elapsed_seconds=time.monotonic() - started, exit_code=0,
                      standard_perambulator_SHA256=projected['output_SHA256'],
                      maximum_full_relative=accepted['maximum_full_relative'],
                      full_gamma5_relative_Frobenius=projected['full_gamma5_relative_Frobenius'],
                      complete_evidence_saved=True, point_intermediate_retired=True)
        write_json(report / 'queue-complete.json', record)
        save()
    except Exception as error:
        code = error.code if isinstance(error, StageError) else 1
        record.update(status='paused' if code == 75 else 'failed', error=str(error),
                      exit_code=code, finished_UTC=utc_now())
        if started is not None:
            record['elapsed_seconds'] = time.monotonic() - started
        save()
        if isinstance(error, StageError):
            raise SystemExit(code)
        raise
    finally:
        if keeper.poll() is None:
            keeper.terminate()
        keeper.wait()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--time', type=int, choices=range(1, 11), required=True)
    main(parser.parse_args().time)
