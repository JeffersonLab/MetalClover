from pathlib import Path
import subprocess,os,time,json,hashlib
root=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU/metal-bsr-prototype')
state=Path('/private/tmp/mg-null32/state.json')
while True:
 s=json.loads(state.read_text())
 if s['status']!='running':break
 time.sleep(5)
if s['status']!='complete':raise SystemExit('Benchmark incomplete; leaving GPU tests stopped')
out=Path('/private/tmp/metalclover-final12-v2');out.mkdir(exist_ok=True)
status={'status':'running','steps':[]}
def save():
 (out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
def run(name,args,cwd=root):
 status['phase']=name;save();start=time.time()
 with (out/(name+'.log')).open('w') as log:
  code=subprocess.call(args,cwd=cwd,env={**os.environ,'VECLIB_MAXIMUM_THREADS':'1','OMP_NUM_THREADS':'1'},stdout=log,stderr=subprocess.STDOUT)
 status['steps'].append({'name':name,'exit':code,'seconds':time.time()-start});save()
 if code:raise RuntimeError(f'{name}: exit {code}')
try:
 flags=['-std=c++17','-O3','-mcpu=apple-m1','-ffp-contract=fast','-Wall','-Wextra','-fobjc-arc','-fblocks','-DACCELERATE_NEW_LAPACK']
 binary='/private/tmp/metalclover-api-build/low_memory_solve'
 run('driver-build',['xcrun','clang++',*flags,'low_memory_solve.mm','-o',binary,'-framework','Metal','-framework','Foundation','-framework','Accelerate'])
 run('analytic-svd',[binary,'--basis-self-test'])
 run('small-real-mg',[binary,'/private/tmp/mg-lowmem12',str(out/'test'),'mg','2','3,3,3,8','8','128','16','4','--basis','near-null-svd','--null-vectors','24','--null-iters','128','--null-restart','8','--null-warm','32','--null-tol','.001','--null-batch','4'],Path('/private/tmp/mg-null32-runtime'))
 run('original-chroma',['/private/tmp/mg-near-null-verify','/private/tmp/mg-null-metadata/correctness12/action.xml','/private/tmp/mg-lowmem12',str(out/'test'),'metal-mg','2',str(out/'native.csv')])
 hashes=[]
 for rhs in range(2):
  new=out/f'test.metal-mg.rhs{rhs}.bin';old=Path(f'/private/tmp/null-batch12.metal-mg.rhs{rhs}.bin')
  a=hashlib.sha256(new.read_bytes()).hexdigest();b=hashlib.sha256(old.read_bytes()).hexdigest()
  hashes.append({'rhs':rhs,'new_sha256':a,'pre_namespace_sha256':b,'bitwise_identical':a==b})
 (out/'namespace-bitwise.json').write_text(json.dumps(hashes,indent=2)+'\n')
 if not all(h['bitwise_identical'] for h in hashes):raise RuntimeError('Namespace/refactor solution bits differ; investigate')
 gcc=Path('/private/tmp/metalclover-api-build/api_contract_gcc')
 while not gcc.exists():
  data=json.loads(Path('/private/tmp/metalclover-namespace-checks/results.json').read_text())
  if any(r['code'] for r in data):raise RuntimeError('Namespace/client build failed')
  time.sleep(3)
 run('gcc-api-contract',[str(gcc),'/private/tmp/mg-lowmem12.bridge'])
 # Rebuild immutable benchmark snapshot with identical flags and output basename.
 snap=Path('/private/tmp/mg-null32-benchmark-source');dest=Path('/private/tmp/mg-null32-repro');dest.mkdir(exist_ok=True)
 run('benchmark-rebuild',['xcrun','clang++',*flags,'metal-bsr-prototype/low_memory_solve.mm','-o',str(dest/'mg-null32-solve'),'-framework','Metal','-framework','Foundation','-framework','Accelerate'],snap)
 original=hashlib.sha256(Path('/private/tmp/mg-null32-solve').read_bytes()).hexdigest()
 rebuilt=hashlib.sha256((dest/'mg-null32-solve').read_bytes()).hexdigest()
 (out/'benchmark-reproduction.json').write_text(json.dumps({'original_sha256':original,'rebuilt_sha256':rebuilt,'binary_identical':original==rebuilt},indent=2)+'\n')
 status['status']='complete'
except BaseException as e:
 status['status']='failed';status['error']=str(e);raise
finally:save()
