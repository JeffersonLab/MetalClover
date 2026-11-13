from pathlib import Path
import subprocess,os,time,json,hashlib
root=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU/metal-bsr-prototype')
checks=Path('/private/tmp/chroma-batch-validation')
while True:
 state=json.loads((checks/'preflight-status.json').read_text())
 if state['status']!='running':break
 time.sleep(3)
if state['status']!='complete':raise RuntimeError('Correctness preflight did not pass')
out=Path('/private/tmp/metalclover-batch-real12');out.mkdir(exist_ok=True)
manifest=json.loads((checks/'benchmark-source-binary.json').read_text())
for name,entry in manifest.items():
 if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=entry['sha256']:
  raise RuntimeError('Timing source/binary changed: '+name)
args=['python3',str(root/'run_guarded.py'),'--log',str(out/'timing.log'),'--resources',str(out/'resources.csv'),'--max-rss-gib','26','--reserve-disk-gib','4','--','/private/tmp/chroma-metal-validation/enabled/mainprogs/tests/t_metal_clover_batch','/Users/admin/Downloads/download/prop_tests/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime','/private/tmp/mg-null-metadata/correctness12/action.xml',str(out/'run'),'96']
status={'status':'running','argv':args,'environment':{'VECLIB_MAXIMUM_THREADS':'1','OMP_NUM_THREADS':'1'},'started':time.time()}
(out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
try:
 code=subprocess.call(args,env={**os.environ,**status['environment']},cwd=out)
 status['exit']=code;status['status']='complete' if code==0 else 'failed'
 if code:raise RuntimeError('Large batch benchmark failed')
finally:
 status['finished']=time.time();(out/'status.json').write_text(json.dumps(status,indent=2)+'\n')
