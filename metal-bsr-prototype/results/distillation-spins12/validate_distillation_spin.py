from pathlib import Path
import json,os,subprocess,time
root=Path('/private/tmp/chroma-distillation-spin'); tests=Path('/private/tmp/chroma-metal-validation/enabled/mainprogs/tests')
env=os.environ.copy();env['SDKROOT']=subprocess.run(['xcrun','--show-sdk-path'],check=True,text=True,capture_output=True).stdout.strip()
env['OMP_NUM_THREADS']='1';env['VECLIB_MAXIMUM_THREADS']='1'
gauge='/Users/admin/Downloads/download/prop_tests/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime'
action='/private/tmp/mg-null-metadata/correctness12/action.xml'
steps=[('regression-rebuild',['make','-o','Makefile','t_metal_clover'],tests),('benchmark-build',['make','-f','Makefile.batch','-o','Makefile','t_metal_clover_batch'],tests),('regression',[str(tests/'t_metal_clover'),gauge,action,'3,3,3,8'],root),('four-spin-smoke',[str(tests/'t_metal_clover_batch'),gauge,action,str(root/'spins1'),'1'],root)]
record={'status':'running','steps':[]}
for name,args,cwd in steps:
 print(name,flush=True);record['stage']=name;started=time.time()
 (root/'validation-status.json').write_text(json.dumps(record,indent=2)+'\n')
 with (root/(name+'.log')).open('w') as log:
  code=subprocess.run(args,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
 record['steps'].append({'name':name,'command':args,'exit':code,'seconds':time.time()-started})
 if code:
  record['status']='failed';(root/'validation-status.json').write_text(json.dumps(record,indent=2)+'\n');raise SystemExit(code)
record['status']='complete';(root/'validation-status.json').write_text(json.dumps(record,indent=2)+'\n')
print('DISTILLATION_SPIN_VALIDATION_COMPLETE',flush=True)
