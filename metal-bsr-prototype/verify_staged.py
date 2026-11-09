#!/usr/bin/env python3
"""Verify saved real-field solutions using stable copies outside FileProvider.

Staging is outside solver timings. Every copied byte is hashed; the original
fields are left intact. This also avoids macOS executable signature-cache
failures after an iCloud placeholder is restored.
"""
import argparse
import hashlib
import json
from math import prod
import os
from pathlib import Path
import shutil
import struct
import subprocess
import tempfile
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('directory', type=Path)
p.add_argument('--bridge', required=True)
p.add_argument('--rhs', type=int, default=2)
p.add_argument('--executable', type=Path, default=Path('build/performance_chroma'))
a = p.parse_args(); d = a.directory
manifest = {'staging_outside_solver_timings': True, 'original_fields_modified': False, 'files': []}
def stage(src, dst, expected=None):
    src = Path(src)
    initial = src.stat()
    if initial.st_flags & 0x40000000:
        print('HYDRATE '+str(src), flush=True)
        subprocess.run(['/usr/bin/brctl','download',str(src.resolve())], check=True)
        start = time.monotonic()
        while src.stat().st_flags & 0x40000000:
            if time.monotonic()-start > 900:
                raise RuntimeError('FileProvider did not restore '+str(src))
            time.sleep(1)
    digest = hashlib.sha256(); count = 0
    print('STAGE '+str(src), flush=True)
    with src.open('rb') as f, dst.open('xb') as out:
        for block in iter(lambda: f.read(8*1024*1024), b''):
            out.write(block); digest.update(block); count += len(block)
    if count != initial.st_size or (expected is not None and count != expected):
        raise RuntimeError(f'Incomplete source: {src}: {count} bytes, expected {expected or initial.st_size}')
    row = {'source':str(src),'bytes':count,'sha256':digest.hexdigest(),
           'hydrated_placeholder':bool(initial.st_flags & 0x40000000)}
    manifest['files'].append(row)
    (d/'verification-staging.json').write_text(json.dumps(manifest,indent=2)+'\n')
    return row

with tempfile.TemporaryDirectory(prefix='mg-chroma-verify-',dir='/private/tmp') as name:
    temp = Path(name)
    stage(a.executable,temp/'chroma-native')
    (temp/'chroma-native').chmod(0o700)
    stage(d/'action.xml',temp/'action.xml')
    bridge = Path(a.bridge+'.bridge')
    row = stage(bridge,temp/'bridge.bridge')
    with (temp/'bridge.bridge').open('rb') as f:
        header = f.read(76)
    if header[:8] != b'CHRMG002':
        raise RuntimeError('Invalid production bridge')
    shape = struct.unpack_from('<4I',header,8)
    volume = prod(shape); vector_bytes = volume*12*16
    if row['bytes'] != 76+volume*(36+12)*16:
        raise RuntimeError('Incomplete production bridge')
    for rhs in range(a.rhs):
        for method in ['chroma-bicgstab','metal-bicgstab','metal-mg']:
            src = d/f'test.{method}.rhs{rhs}.bin'
            row = stage(src,temp/src.name,vector_bytes)
            old = d/'first-cpu-reference.json'
            if rhs==0 and method=='chroma-bicgstab' and old.exists():
                if row['sha256'] != json.loads(old.read_text())['solution_sha256']:
                    raise RuntimeError('Restored CPU reference differs from its original checkpoint hash')
    log = d/'verify.log'
    if log.exists():
        if 'PASS mode=verify' in log.read_text():
            raise RuntimeError('Native verification already passed; refusing to repeat it')
        n = 1
        while (d/f'verify-attempt{n}.log').exists(): n += 1
        log.rename(d/f'verify-attempt{n}.log')
    cmd = [str(temp/'chroma-native'),str(temp/'action.xml'),str(temp/'bridge'),str(temp/'test'),'verify',str(a.rhs)]
    manifest['native_argv'] = cmd
    print('VERIFY '+repr(cmd),flush=True)
    with log.open('x') as out:
        result = subprocess.run(cmd,env=dict(os.environ,OMP_NUM_THREADS='8'),stdout=out,stderr=subprocess.STDOUT)
    manifest['native_returncode'] = result.returncode
    (d/'verification-staging.json').write_text(json.dumps(manifest,indent=2)+'\n')
    if result.returncode:
        raise RuntimeError('Native verification failed; see '+str(log))
    shutil.copyfile(temp/'test.chroma-checks.csv',d/'test.chroma-checks.csv')
    with (d/'files.sha256').open('w') as out:
        for row in manifest['files']:
            out.write(row['sha256']+'  '+row['source']+'\n')
print('PASS staged original Chroma FP64 verification; all copies hashed',flush=True)
