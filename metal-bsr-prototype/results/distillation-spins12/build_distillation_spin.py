from pathlib import Path
import datetime,json,os,subprocess,time
out=Path('/private/tmp/chroma-distillation-spin')
build=Path('/private/tmp/chroma-metal-validation/enabled')
env=os.environ.copy()
env['SDKROOT']=subprocess.run(['xcrun','--show-sdk-path'],check=True,capture_output=True,text=True).stdout.strip()
env['OMP_NUM_THREADS']='1';env['VECLIB_MAXIMUM_THREADS']='1'
steps=[('native-objects',['make','-o','Makefile','-j1','util/ferm/mgproton.o','meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.o'],build/'lib'),
       ('native-archive',['/usr/bin/ar','r','libchroma.a','util/ferm/mgproton.o','meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.o'],build/'lib'),
       ('native-index',['/usr/bin/ranlib','libchroma.a'],build/'lib'),
       ('regression-build',['make','-o','Makefile','t_metal_clover'],build/'mainprogs/tests'),
       ('benchmark-build',['make','-f','Makefile.batch','-o','Makefile','t_metal_clover_batch'],build/'mainprogs/tests')]
records=[]
for name,cmd,cwd in steps:
 print(name,flush=True);started=time.time()
 with (out/(name+'.log')).open('w') as log:
  code=subprocess.run(cmd,cwd=cwd,env=env,stdout=log,stderr=subprocess.STDOUT).returncode
 records.append({'name':name,'command':cmd,'cwd':str(cwd),'exit':code,'seconds':time.time()-started})
 (out/'build-status.json').write_text(json.dumps(records,indent=2)+'\n')
 if code: raise SystemExit(code)
print('DISTILLATION_BUILD_COMPLETE',flush=True)
