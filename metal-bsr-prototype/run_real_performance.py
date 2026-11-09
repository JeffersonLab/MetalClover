#!/usr/bin/env python3
"""Sequential real-field native Chroma/Metal timing and FP64 verification.

The action XML, input.json and native production bridge must already exist.
--wait-pid/--skip-cpu can continue a CPU run started with run_guarded.py.
"""
import argparse
import atexit
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import shutil
import tempfile
import time

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('directory', type=Path)
p.add_argument('--bridge', required=True)
p.add_argument('--rhs', type=int, default=2)
p.add_argument('--blocks', default='4,4,4,8')
p.add_argument('--skip-cpu', action='store_true')
p.add_argument('--wait-pid', type=int)
p.add_argument('--chroma-config', default='/Users/admin/Downloads/macos-chroma/install/chroma/bin/chroma-config')
p.add_argument('--chroma-source', default='/Users/admin/Downloads/qcd-macos-share-20261005/chroma.devel')
a = p.parse_args(); d = a.directory; prefix = str(d/'test')
info = json.loads((d/'input.json').read_text())
commands = []
def run(command, logfile=None, env=None):
    print('RUN '+repr(command), flush=True)
    commands.append({'argv': list(map(str,command)), 'log': str(logfile) if logfile else None, 'env': env or {}})
    (d/'commands.json').write_text(json.dumps(commands, indent=2)+'\n')
    active_env = dict(os.environ, **(env or {}))
    if logfile:
        with Path(logfile).open('x') as f:
            subprocess.run(list(map(str,command)), env=active_env, stdout=f, stderr=subprocess.STDOUT, check=True)
    else:
        subprocess.run(list(map(str,command)), env=active_env, check=True)
if a.wait_pid:
    if not a.skip_cpu:
        p.error('--wait-pid requires --skip-cpu')
    print(f'Waiting for CPU PID {a.wait_pid}', flush=True)
    while True:
        try:
            os.kill(a.wait_pid, 0)
        except ProcessLookupError:
            break
        time.sleep(5)
run(['make','-f','Makefile.production','build/performance_compare','build/performance_chroma','build/correctness_chroma',
     'CHROMA_CONFIG='+a.chroma_config,'CHROMA_SOURCE='+a.chroma_source], d/'benchmark-build.log')
# iCloud can dehydrate a build product while this long benchmark is running.
# Fresh, byte-identical executable copies avoid the macOS signature-cache
# rejection observed when launching a restored FileProvider executable.
images = tempfile.TemporaryDirectory(prefix='mg-performance-exec-',dir='/private/tmp')
atexit.register(images.cleanup)
executables = {}; image_records = []
for name in ['performance_compare','performance_chroma','correctness_chroma']:
    src = Path('build')/name
    if src.stat().st_flags & 0x40000000:
        subprocess.run(['/usr/bin/brctl','download',str(src.resolve())],check=True)
        start = time.monotonic()
        while src.stat().st_flags & 0x40000000:
            if time.monotonic()-start>900: raise RuntimeError('Cannot restore '+str(src))
            time.sleep(1)
    dst = Path(images.name)/name
    shutil.copyfile(src,dst);dst.chmod(0o700)
    if dst.stat().st_size!=src.stat().st_size:
        raise RuntimeError('Incomplete executable copy: '+str(src))
    source_hash = hashlib.sha256(src.read_bytes()).hexdigest()
    assert source_hash == hashlib.sha256(dst.read_bytes()).hexdigest()
    executables[name] = str(dst)
    image_records.append({'build_product':str(src),'executed_copy':str(dst),'sha256':source_hash})
(d/'executables.json').write_text(json.dumps(image_records,indent=2)+'\n')
if not a.skip_cpu:
    run([sys.executable,'run_guarded.py','--log',str(d/'chroma.log'),'--resources',str(d/'chroma-resources.csv'),'--',
         executables['performance_chroma'],str(d/'action.xml'),a.bridge,prefix,'cpu',str(a.rhs)], env={'OMP_NUM_THREADS':'8'})
with (d/'test.chroma-solves.csv').open() as f:
    solved = list(csv.DictReader(f))
assert len(solved)==a.rhs and all(float(row['full_relative'])<=1e-12 for row in solved)
assert 'PASS mode=cpu' in (d/'chroma.log').read_text()
run([executables['correctness_chroma'],info['gauge'],str(d/'action.xml'),a.bridge,prefix,'state-only'],
    d/'stout-state-check.log', {'OMP_NUM_THREADS':'8'})
run([sys.executable,'run_guarded.py','--log',str(d/'metal.log'),'--resources',str(d/'metal-resources.csv'),'--',
     executables['performance_compare'],a.bridge,prefix,str(a.rhs),a.blocks], env={'VECLIB_MAXIMUM_THREADS':'1'})
run([sys.executable,'verify_staged.py',str(d),'--bridge',a.bridge,'--rhs',str(a.rhs)],d/'verification-stage.log')
run([sys.executable,'performance_report.py',str(d),'--rhs',str(a.rhs)], d/'report.log')
print('COMPLETE all timings and independent full native FP64 checks', flush=True)
