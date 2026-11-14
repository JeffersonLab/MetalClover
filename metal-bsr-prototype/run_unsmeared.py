#!/usr/bin/env python3
"""Run/check a checkpointed native Chroma unsmeared propagator job."""
import argparse
import ast
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import signal
import struct
import subprocess
import sys
import time
import xml.etree.ElementTree as ET


def atomic_json(path, data):
    temp = Path(str(path) + '.tmp')
    temp.write_text(json.dumps(data, indent=2) + '\n')
    temp.replace(path)


def sha256(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def checkpoint(prefix):
    path = Path(str(prefix) + '.state')
    if not path.exists():
        return {}
    rows = path.read_text().splitlines()
    if len(rows) != 8 or rows[0] != 'METAL_UNSMEARED_V1':
        raise RuntimeError('invalid checkpoint')
    return {'completed_rhs': int(rows[2]), 'max_full_relative': float(rows[4]),
            'qprop_seconds': float(rows[5]), 'verify_seconds': float(rows[6]),
            'write_seconds': float(rows[7])}


def read_npy_header(path):
    with Path(path).open('rb') as f:
        if f.read(8) != b'\x93NUMPY\x01\x00':
            raise RuntimeError('not a version-1 NPY file')
        length = struct.unpack('<H', f.read(2))[0]
        header = ast.literal_eval(f.read(length).decode('ascii').strip())
    return header, 10 + length


def summarize(prefix):
    meta = ET.parse(str(prefix) + '.meta.xml').getroot()
    count = int(meta.findtext('num_vectors'))
    first = int(meta.findtext('first_vector'))
    dims = tuple(map(int, meta.findtext('lattSize').split()))
    shape = (count, 4, dims[2], dims[1], dims[0], 3, 4)
    state = checkpoint(prefix)
    if state['completed_rhs'] != 4 * count:
        raise RuntimeError('output is incomplete')
    header, offset = read_npy_header(str(prefix) + '.npy')
    if header != {'descr': '<c16', 'fortran_order': False, 'shape': shape}:
        raise RuntimeError('NPY dtype/shape mismatch')
    expected_bytes = offset + math.prod(shape) * 16
    if Path(str(prefix) + '.npy').stat().st_size != expected_bytes:
        raise RuntimeError('NPY byte size mismatch')
    with Path(str(prefix) + '.columns.csv').open() as f:
        columns = list(csv.DictReader(f))
    if len(columns) != count * 4:
        raise RuntimeError('column journal length mismatch')
    for i, row in enumerate(columns):
        vector, spin = i % count, i // count
        if (int(row['ordered_rhs']), int(row['vector']), int(row['source_spin']), int(row['column'])) != (i, first + vector, spin, 4 * vector + spin):
            raise RuntimeError('journal vector/spin ordering mismatch')
        for key in ['full_relative', 'qprop_relative']:
            if not 0 <= float(row[key]) <= 1e-12:
                raise RuntimeError('full FP64 residual acceptance failed')
        if not (math.isfinite(float(row['source_norm'])) and float(row['source_norm']) > 0):
            raise RuntimeError('invalid eigenvector source norm')
    sessions = sorted(prefix.parent.glob(prefix.name + '.session-*.status.json'))
    wall = sum(json.loads(p.read_text()).get('elapsed_seconds', 0) for p in sessions)
    flops = 0
    hierarchies = 0
    for path in sorted(prefix.parent.glob(prefix.name + '.session-*.log')):
        for line in path.read_text().splitlines():
            if ' array: rhs=' in line:
                m = re.search(r'nominal_operator_flops=(\d+)', line)
                if m:
                    flops += int(m.group(1))
            if line.startswith('NEAR_NULL_SVD '):
                hierarchies += 1
    result = {'status': 'complete', 'method': meta.findtext('method'), 'source_time': 0,
              'sink_time': 0, 'lattice': dims, 'first_vector': first, 'num_vectors': count,
              'source_spins': 4, 'physical_rhs': 4 * count, 'numpy_shape': shape,
              'numpy_axes': meta.findtext('axes_C_order'), 'dtype': '<c16',
              'input_eigen_precision_bits': int(meta.findtext('input_eigen_precision_bits')),
              'full_CPU_FP64_tolerance': 1e-12,
              'maximum_full_relative': max(float(x['full_relative']) for x in columns),
              'maximum_qprop_relative': max(float(x['qprop_relative']) for x in columns),
              'iterations_min': min(int(x['iterations']) for x in columns),
              'iterations_max': max(int(x['iterations']) for x in columns),
              'iterations_sum': sum(int(x['iterations']) for x in columns),
              'MG_hierarchy_builds_across_sessions': hierarchies,
              'session_count': len(sessions), 'elapsed_seconds_across_sessions': wall,
              **state, 'payload_bytes': expected_bytes,
              'payload_sha256': sha256(str(prefix) + '.npy'),
              'nominal_dispatched_operator_flops': flops,
              'nominal_operator_GFLOP_per_second_total_wall': flops / (1e9 * wall) if wall else None,
              'rate_scope': 'nominal dispatched GPU operator work, including inactive/padded lanes; CPU FP64 and Krylov-vector arithmetic excluded from numerator; session wall includes setup, input, residual checks and output, excluding final checksum reread and SHA256',
              'source_order': meta.findtext('source_order'),
              'completed_payload_checksums': 'all column FNV1a64 checks independently reread by native --resume verifier'}
    atomic_json(Path(str(prefix) + '.complete.json'), result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--gauge', type=Path, required=True)
    parser.add_argument('--action', type=Path, required=True)
    parser.add_argument('--eigenvectors', type=Path, required=True)
    parser.add_argument('--prefix', type=Path, required=True)
    parser.add_argument('--first', type=int, default=0)
    parser.add_argument('--vectors', type=int, default=5184)
    parser.add_argument('--max-rhs', type=int, default=16)
    parser.add_argument('--method', choices=['MetalMGClover', 'MetalBiCGClover'], default='MetalMGClover')
    parser.add_argument('--resume', action='store_true')
    parser.add_argument('--max-rss-gib', type=float, default=26)
    parser.add_argument('--reserve-disk-gib', type=float, default=4)
    a = parser.parse_args()
    a.prefix = a.prefix.resolve()
    a.prefix.parent.mkdir(parents=True, exist_ok=True)
    prior = list(a.prefix.parent.glob(a.prefix.name + '.session-*.status.json'))
    index = len(prior)
    session = Path(str(a.prefix) + f'.session-{index:03d}')
    status_path = Path(str(a.prefix) + '.run.json')
    command = [str(a.binary.resolve()), str(a.gauge.resolve()), str(a.action.resolve()),
               str(a.eigenvectors.resolve()), str(a.prefix), str(a.first), str(a.vectors),
               str(a.max_rhs), a.method]
    if a.resume:
        command.append('--resume')
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    record = {'status': 'running', 'command': command, 'total_rhs': a.vectors * 4,
              'started_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
              'binary_sha256': sha256(a.binary), 'inputs_sha256': {str(p.resolve()): sha256(p) for p in [a.gauge, a.action, a.eigenvectors]},
              'environment': {'OMP_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '1'},
              'limits': {'max_rss_gib': a.max_rss_gib, 'reserve_disk_gib': a.reserve_disk_gib}}
    requested = []
    def on_signal(sig, _):
        requested.append(sig)
    signal.signal(signal.SIGTERM, on_signal)
    signal.signal(signal.SIGINT, on_signal)
    started = time.monotonic()
    reason = ''
    stop_started = None
    peak_rss = 0.0
    minimum_free = float('inf')
    with Path(str(session) + '.log').open('x') as log, Path(str(session) + '.resources.csv').open('x', newline='') as resources:
        writer = csv.writer(resources)
        writer.writerow(['elapsed_seconds', 'rss_gib', 'disk_free_gib', 'completed_rhs'])
        child = subprocess.Popen(command, cwd=a.prefix.parent, env=env, stdout=log, stderr=subprocess.STDOUT)
        keeper = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(child.pid)])
        record['PID'] = child.pid
        atomic_json(status_path, record)
        while child.poll() is None:
            elapsed = time.monotonic() - started
            proc = subprocess.run(['ps', '-o', 'rss=', '-p', str(child.pid)], text=True, capture_output=True)
            if proc.returncode and child.poll() is None:
                reason = 'process RSS monitoring failed: ' + proc.stderr.strip()
            rss = int(proc.stdout.strip() or '0') / 1024**2
            free = shutil.disk_usage(a.prefix.parent).free / 1024**3
            peak_rss = max(peak_rss, rss)
            minimum_free = min(minimum_free, free)
            current = checkpoint(a.prefix)
            writer.writerow([elapsed, rss, free, current.get('completed_rhs', 0)])
            resources.flush()
            if not reason and (rss > a.max_rss_gib or free < a.reserve_disk_gib):
                reason = f'resource guard: RSS={rss:.3f} GiB free_disk={free:.3f} GiB'
            if requested and not reason:
                reason = 'requested pause'
            if reason and stop_started is None:
                child.send_signal(signal.SIGTERM)
                stop_started = time.monotonic()
            if stop_started is not None and time.monotonic() - stop_started > 90:
                child.kill()
            record.update(current)
            record.update(elapsed_seconds=elapsed, RSS_GiB=rss, free_disk_GiB=free,
                          peak_RSS_GiB=peak_rss, minimum_free_disk_GiB=minimum_free)
            if reason:
                record['pause_reason'] = reason
            atomic_json(status_path, record)
            time.sleep(5)
        code = child.wait()
        try:
            keeper.wait(timeout=5)
        except subprocess.TimeoutExpired:
            keeper.terminate()
            keeper.wait()
    record.update(checkpoint(a.prefix))
    record.update(exit_code=code, elapsed_seconds=time.monotonic() - started,
                  peak_RSS_GiB=peak_rss, minimum_free_disk_GiB=minimum_free,
                  finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
    record['status'] = 'complete' if code == 0 and record.get('completed_rhs') == a.vectors * 4 else ('paused' if code == 75 else 'failed')
    if reason:
        record['pause_reason'] = reason
    atomic_json(Path(str(session) + '.status.json'), record)
    atomic_json(status_path, record)
    if record['status'] == 'complete':
        record['status'] = 'verifying'
        atomic_json(status_path, record)
        verify_command = command if '--resume' in command else command + ['--resume']
        with Path(str(session) + '.verify.log').open('x') as log:
            verify = subprocess.run(verify_command, cwd=a.prefix.parent, env=env, stdout=log, stderr=subprocess.STDOUT)
        if verify.returncode:
            record['status'] = 'verification_failed'
            record['verify_exit'] = verify.returncode
            atomic_json(status_path, record)
            return 1
        try:
            result = summarize(a.prefix)
        except Exception as error:
            record['status'] = 'verification_failed'
            record['verify_error'] = str(error)
            atomic_json(status_path, record)
            raise
        record['status'] = 'complete'
        record['payload_sha256'] = result['payload_sha256']
        atomic_json(status_path, record)
        print(json.dumps(result, indent=2), flush=True)
        return 0
    print(json.dumps(record, indent=2), flush=True)
    return code or 1


if __name__ == '__main__':
    sys.exit(main())
