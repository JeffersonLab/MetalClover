#!/usr/bin/env python3
"""Finish the sequential real 32³ comparison after an already running BiCGStab.

Native verification and GPU inversions never overlap. The controller owns
later process groups so an explicit pause can stop all of their children.
It refuses to accept incomplete solves or native residuals above 1e-12.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import time

PROJECT=Path.cwd()
D=(PROJECT/'results/performance32-fit')
RUNTIME=Path('/private/tmp/mg-fit32-runtime')
RUNTIME.mkdir(exist_ok=True)
PREFIX='/private/tmp/mg-fit32-test'
EXE='/private/tmp/mg-lowmem-solve-v3'
VERIFY='/private/tmp/mg-lowmem-verify'
# Pin small runtime inputs outside FileProvider; the large bridge/fields are
# already local scratch files. This avoids stalled reads of cloud placeholders.
for src,dest in [*[(p,RUNTIME/p.name) for p in PROJECT.glob('*.metal')],(PROJECT/'run_guarded.py',RUNTIME/'run_guarded.py'),(D/'action.xml',RUNTIME/'action32.xml'),(PROJECT/'results/correctness12/action.xml',RUNTIME/'action12.xml')]:
    if src.stat().st_flags & 0x40000000:
        subprocess.run(['brctl','download',str(src)],check=True)
        start=time.monotonic()
        while src.stat().st_flags & 0x40000000:
            if time.monotonic()-start>900:raise RuntimeError('FileProvider could not restore '+str(src))
            time.sleep(1)
    shutil.copyfile(src,dest)
with (D/'runtime-files.sha256').open('w') as out:
    for src in sorted(RUNTIME.iterdir()):
        if src.is_file():out.write(hashlib.sha256(src.read_bytes()).hexdigest()+'  '+str(src)+'\n')
child=None
state={'phase':'waiting for existing BiCGStab','pid':os.getpid()}

def save(**values):
    state.update(values);(D/'controller-state.json').write_text(json.dumps(state,indent=2)+'\n')

def stop(signum,frame):
    if child is not None and child.poll() is None:
        os.killpg(child.pid,signal.SIGTERM)
    save(phase='stopped by signal',signal=signum)
    raise SystemExit(128+signum)

signal.signal(signal.SIGTERM,stop);signal.signal(signal.SIGINT,stop)

def rows(path):
    with Path(path).open(newline='') as f:return list(csv.DictReader(f))

def accept(path,method,count):
    data=rows(path)
    if len(data)!=count or sorted(int(r['rhs']) for r in data)!=list(range(count)):
        raise RuntimeError('missing complete fields in '+str(path))
    for r in data:
        if r['method']!=method or not 0<=float(r['full_relative'])<=1e-12:
            raise RuntimeError('full FP64 residual rejected: '+repr(r))
    return data

def run(name,args,env=None,rss=26):
    global child
    save(phase=name,argv=args)
    command=['python3',str(RUNTIME/'run_guarded.py'),'--log',str(D/(name+'.log')),'--resources',str(D/(name+'-resources.csv')),'--max-rss-gib',str(rss),'--reserve-disk-gib','2','--',*args]
    child=subprocess.Popen(command,env=dict(os.environ,**(env or {})),start_new_session=True,cwd=RUNTIME)
    save(child_group=child.pid)
    code=child.wait();child=None
    if code:raise RuntimeError(f'{name} failed ({code})')

save()
while 'PASS all full FP64 residuals' not in (D/'bicg.log').read_text():
    listing=subprocess.check_output(['ps','-ax','-o','comm='],text=True)
    if '/private/tmp/mg-lowmem-solve-v2' not in listing:
        raise RuntimeError('existing BiCGStab exited without a complete pass')
    time.sleep(5)
accept(PREFIX+'.bicg.csv','metal-bicgstab',2)
for ext in ['csv','trace.csv']:shutil.copyfile(PREFIX+'.bicg.'+ext,D/('test.bicg.'+ext))
run('verify-bicg',[VERIFY,str(RUNTIME/'action32.xml'),'/private/tmp/mg-fit32',PREFIX,'metal-bicgstab','2',str(D/'test.bicg-native.csv')],{'OMP_NUM_THREADS':'8'})
accept(D/'test.bicg-native.csv','metal-bicgstab',2)
# Validate the stronger coarse basis / bottom cap through the complete small
# real-field solve before applying that preconditioner to the large field.
SMALL='/private/tmp/mg-lowmem12-v3'
run('small-mg-v3',[EXE,'/private/tmp/mg-lowmem12',SMALL,'mg','1','3,3,3,8','8','256','64','8'],{'VECLIB_MAXIMUM_THREADS':'1'})
accept(SMALL+'.mg.csv','metal-mg',1)
run('verify-small-mg-v3',[VERIFY,str(RUNTIME/'action12.xml'),'/private/tmp/mg-lowmem12',SMALL,'metal-mg','1',str(D/'small-mg-v3-native.csv')],{'OMP_NUM_THREADS':'8'})
accept(D/'small-mg-v3-native.csv','metal-mg',1)
run('mg',[EXE,'/private/tmp/mg-fit32',PREFIX,'mg','2','4,4,4,4','8','256','64','8'],{'VECLIB_MAXIMUM_THREADS':'1'})
accept(PREFIX+'.mg.csv','metal-mg',2)
for ext in ['csv','trace.csv']:shutil.copyfile(PREFIX+'.mg.'+ext,D/('test.mg.'+ext))
run('verify-mg',[VERIFY,str(RUNTIME/'action32.xml'),'/private/tmp/mg-fit32',PREFIX,'metal-mg','2',str(D/'test.mg-native.csv')],{'OMP_NUM_THREADS':'8'})
accept(D/'test.mg-native.csv','metal-mg',2)
save(phase='hashing accepted fields')
files=[]
for src in [Path('/private/tmp/mg-fit32.bridge'),*[Path(PREFIX+f'.{method}.rhs{rhs}.bin') for method in ['metal-bicgstab','metal-mg'] for rhs in range(2)]]:
    h=hashlib.sha256();n=0
    with src.open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b);n+=len(b)
    files.append({'path':str(src),'bytes':n,'sha256':h.hexdigest()})
(D/'accepted-files.json').write_text(json.dumps(files,indent=2)+'\n')
with (D/'accepted-files.sha256').open('w') as f:
    for r in files:f.write(r['sha256']+'  '+r['path']+'\n')
save(phase='complete',accepted_full_tolerance=1e-12)
print('PASS complete original-native FP64 verified 32^3x256 comparison',flush=True)
