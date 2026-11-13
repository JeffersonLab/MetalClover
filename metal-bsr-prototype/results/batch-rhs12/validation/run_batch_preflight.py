from pathlib import Path
import subprocess,os,time,json
out=Path('/private/tmp/chroma-batch-validation')
root=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU/metal-bsr-prototype')
env={**os.environ,'VECLIB_MAXIMUM_THREADS':'1','OMP_NUM_THREADS':'1'}
state={'status':'running','steps':[]}
def save(): (out/'preflight-status.json').write_text(json.dumps(state,indent=2)+'\n')
def run(name,args):
 state['phase']=name;save();start=time.time()
 with (out/(name+'.log')).open('w') as log:
  code=subprocess.call(args,env=env,stdout=log,stderr=subprocess.STDOUT)
 state['steps'].append({'name':name,'exit':code,'seconds':time.time()-start});save()
 print(name,code,flush=True)
 if code:raise RuntimeError(name+' failed')
try:
 run('real12-regression',['/private/tmp/chroma-metal-validation/enabled/mainprogs/tests/t_metal_clover','/Users/admin/Downloads/download/prop_tests/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime','/private/tmp/mg-null-metadata/correctness12/action.xml','3,3,3,8'])
 run('gcc-contract',['/private/tmp/metalclover-batch-build/api_contract_gcc','/private/tmp/mg-lowmem12.bridge','bicg'])
 binary=Path('/private/tmp/chroma-metal-validation/enabled/mainprogs/tests/t_metal_clover_batch')
 while not binary.exists():time.sleep(2)
 run('benchmark-smoke',[str(binary),'/Users/admin/Downloads/download/prop_tests/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime','/private/tmp/mg-null-metadata/correctness12/action.xml',str(out/'smoke'),'2'])
 state['status']='complete'
except BaseException as e:state['status']='failed';state['error']=str(e);raise
finally:save()
