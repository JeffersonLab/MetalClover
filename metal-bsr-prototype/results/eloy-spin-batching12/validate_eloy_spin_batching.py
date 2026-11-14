from pathlib import Path
import json,os,subprocess,time
root=Path('/private/tmp/chroma-eloy-spin-batching');tests=Path('/private/tmp/chroma-metal-validation/enabled/mainprogs/tests')
env=os.environ.copy();env['SDKROOT']=subprocess.run(['xcrun','--show-sdk-path'],check=True,text=True,capture_output=True).stdout.strip()
env['OMP_NUM_THREADS']='1';env['VECLIB_MAXIMUM_THREADS']='1'
gauge='/Users/admin/Downloads/download/prop_tests/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime';action='/private/tmp/mg-null-metadata/correctness12/action.xml'
steps=[('regression-build',['make','-o','Makefile','t_metal_clover'],tests),('benchmark-build',['make','-f','Makefile.batch','-o','Makefile','t_metal_clover_batch'],tests),('regression',[str(tests/'t_metal_clover'),gauge,action,'3,3,3,8'],root),('spin-major-smoke',[str(tests/'t_metal_clover_batch'),gauge,action,str(root/'spins1'),'1'],root)]
record={'status':'running','steps':[]}
for name,cmd,cwd in steps:
 print(name,flush=True);record['stage']=name;started=time.time()
 (root/'validation-status.json').write_text(json.dumps(record,indent=2)+'\n')
 with (root/(name+'.log')).open('w') as log:
  result=subprocess.run(cmd,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT)
 record['steps'].append({'name':name,'command':cmd,'exit':result.returncode,'seconds':time.time()-started})
 if result.returncode:
  record['status']='failed';(root/'validation-status.json').write_text(json.dumps(record,indent=2)+'\n');raise SystemExit(result.returncode)
record['status']='complete';(root/'validation-status.json').write_text(json.dumps(record,indent=2)+'\n')
print('ELOY_SPIN_BATCH_VALIDATION_COMPLETE',flush=True)
