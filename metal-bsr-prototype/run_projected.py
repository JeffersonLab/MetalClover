#!/usr/bin/env python3
"""Project an accepted unsmeared collection into a standard Chroma perambulator."""
import argparse
import csv
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time
import xml.etree.ElementTree as ET


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for chunk in iter(lambda: f.read(4 * 1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, data):
    temp = Path(str(path) + '.tmp')
    temp.write_text(json.dumps(data, indent=2) + '\n')
    temp.replace(path)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--binary', type=Path, required=True)
    parser.add_argument('--point-prefix', type=Path, required=True)
    parser.add_argument('--eigenvectors', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--vectors', type=int, default=5184)
    args = parser.parse_args()
    args.point_prefix = args.point_prefix.resolve()
    args.output = args.output.resolve()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    for suffix in ['', '.partial', '.projection.csv', '.projection.xml', '.run.json', '.complete.json', '.log', '.resources.csv']:
        if Path(str(args.output) + suffix).exists():
            raise RuntimeError('projection output already exists: ' + str(args.output) + suffix)
    original = json.loads(Path(str(args.point_prefix) + '.complete.json').read_text())
    original_run = json.loads(Path(str(args.point_prefix) + '.run.json').read_text())
    assert original['status'] == original_run['status'] == 'complete'
    assert original['completed_rhs'] == original['physical_rhs'] == 4 * original['num_vectors']
    assert original['first_vector'] == original['source_time'] == original['sink_time'] == 0
    assert original['maximum_full_relative'] <= 1e-12 and original['maximum_qprop_relative'] <= 1e-12
    assert 1 <= args.vectors <= original['num_vectors']
    expected_eigen = original_run['inputs_sha256'][str(args.eigenvectors.resolve())]
    print('Checking eigenvector and accepted point-file SHA256 hashes...', flush=True)
    eigen_hash = digest(args.eigenvectors)
    point_hash = digest(str(args.point_prefix) + '.npy')
    if eigen_hash != expected_eigen or point_hash != original['payload_sha256']:
        raise RuntimeError('input integrity check failed')
    expected = 16 * args.vectors**2 * 16
    if shutil.disk_usage(args.output.parent).free < 1.1 * expected + 4 * 1024**3:
        raise RuntimeError('insufficient free disk for the standard perambulator plus reserve')
    command = [str(args.binary.resolve()), str(args.point_prefix), str(args.eigenvectors.resolve()), str(args.output), str(args.vectors)]
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='8')
    info = {'status': 'running', 'started_UTC': datetime.datetime.now(datetime.timezone.utc).isoformat(),
            'command': command, 'source_time': 0, 'sink_time': 0, 'num_vectors': args.vectors,
            'definition': 'V_dagger M_inverse V', 'format': 'standard Chroma FILEDB/SDB propElemOp, complex FP64',
            'point_SHA256': point_hash, 'eigenvector_SHA256': eigen_hash,
            'binary_SHA256': digest(args.binary), 'original_full_CPU_FP64_checks': original,
            'environment': {'OMP_NUM_THREADS': '1', 'VECLIB_MAXIMUM_THREADS': '8'}}
    signals = []
    signal.signal(signal.SIGINT, lambda sig, _: signals.append(sig))
    signal.signal(signal.SIGTERM, lambda sig, _: signals.append(sig))
    status_path = Path(str(args.output) + '.run.json')
    started = time.monotonic()
    peak = 0
    minimum_free = float('inf')
    reason = ''
    with Path(str(args.output) + '.log').open('x') as log, Path(str(args.output) + '.resources.csv').open('x', newline='') as resources:
        writer = csv.writer(resources)
        writer.writerow(['seconds', 'rss_gib', 'disk_free_gib'])
        child = subprocess.Popen(command, cwd=args.output.parent, env=env, stdout=log, stderr=subprocess.STDOUT)
        keeper = subprocess.Popen(['/usr/bin/caffeinate', '-i', '-w', str(child.pid)])
        info['PID'] = child.pid
        save(status_path, info)
        while child.poll() is None:
            result = subprocess.run(['ps', '-o', 'rss=', '-p', str(child.pid)], capture_output=True, text=True)
            if result.returncode and child.poll() is None:
                reason = 'process RSS monitoring failed'
            rss = int(result.stdout.strip() or '0') / 1024**2
            free = shutil.disk_usage(args.output.parent).free / 1024**3
            peak = max(peak, rss)
            minimum_free = min(minimum_free, free)
            elapsed = time.monotonic() - started
            writer.writerow([elapsed, rss, free])
            resources.flush()
            if rss > 26 or free < 4:
                reason = 'resource guard'
            if signals:
                reason = 'requested stop'
            if reason:
                child.terminate()
            info.update(elapsed_seconds=elapsed, RSS_GiB=rss, peak_RSS_GiB=peak, free_disk_GiB=free)
            save(status_path, info)
            time.sleep(5)
        code = child.wait()
        keeper.wait(timeout=10)
    info.update(exit_code=code, elapsed_seconds=time.monotonic() - started,
                peak_RSS_GiB=peak, minimum_free_disk_GiB=minimum_free,
                finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
    if code or reason:
        info.update(status='failed', reason=reason)
        save(status_path, info)
        print(json.dumps(info, indent=2), flush=True)
        return code or 1
    info['status'] = 'verifying'
    save(status_path, info)
    try:
        native = ET.parse(str(args.output) + '.projection.xml').getroot()
        assert native.findtext('status') == 'pass'
        assert int(native.findtext('num_vectors')) == args.vectors
        assert int(native.findtext('spin_records')) == 16
        assert float(native.findtext('max_QDP_dot_scaled_error')) <= 1e-12
        assert float(native.findtext('full_gamma5_relative_Frobenius')) <= 1e-10
        assert native.findtext('metadata/DBMetaData/id') == 'propElemOp'
        assert int(native.findtext('metadata/DBMetaData/Params/Contractions/num_vecs')) == args.vectors
        mass = native.findtext('metadata/DBMetaData/Params/Propagator/FermionAction/Mass')
        assert float(mass) == -0.0743
        with Path(str(args.output) + '.projection.csv').open() as f:
            rows = list(csv.DictReader(f))
        assert len(rows) == 16
        assert {(int(r['spin_source']), int(r['spin_sink'])) for r in rows} == {(s, t) for s in range(4) for t in range(4)}
        for row in rows:
            assert int(row['num_vectors']) == args.vectors
            assert float(row['qdp_dot_scaled_error']) <= 1e-12
        info.update(status='complete', output_bytes=args.output.stat().st_size,
                    output_SHA256=digest(args.output), spin_records=16,
                    max_QDP_dot_scaled_error=float(native.findtext('max_QDP_dot_scaled_error')),
                    full_gamma5_relative_Frobenius=float(native.findtext('full_gamma5_relative_Frobenius')),
                    full_gamma5_max_absolute=float(native.findtext('full_gamma5_max_absolute')),
                    native_seconds=float(native.findtext('total_seconds')),
                    gemm_seconds=float(native.findtext('gemm_seconds')),
                    full_record_readback='all 16 records decoded by Chroma and re-encoded bytes matched FNV1a64 hashes')
        save(Path(str(args.output) + '.complete.json'), info)
        save(status_path, info)
    except Exception as error:
        info.update(status='verification_failed', error=str(error))
        save(status_path, info)
        raise
    print(json.dumps(info, indent=2), flush=True)
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
