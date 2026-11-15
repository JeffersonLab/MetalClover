from pathlib import Path
import subprocess,os,json,signal,time,datetime,hashlib
root=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU');scratch=Path('/private/tmp/chroma-time1-production')
bin=Path('/private/tmp/chroma-metal-validation/enabled/mainprogs/tests')
point=Path('/Users/admin/Downloads/MetalClover-results/cfg1000-t1-unsmeared')
peram=Path('/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_1.sdb1000')
gauge=Path('/Users/admin/Downloads/download/prop_tests/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime')
eigen=Path('/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_1.mod1000')
action=root/'metal-bsr-prototype/results/perambulator-t1-5184/action.xml'
commands=[('solve',['python3',root/'metal-bsr-prototype/run_unsmeared.py','--binary',bin/'t_metal_clover_unsmeared_time','--gauge',gauge,'--action',action,'--eigenvectors',eigen,'--prefix',point,'--first','0','--vectors','5184','--max-rhs','16','--method','MetalMGClover','--time','1']),('project',['python3',root/'metal-bsr-prototype/run_projected.py','--binary',bin/'project_perambulator_time','--point-prefix',point,'--eigenvectors',eigen,'--output',peram,'--vectors','5184']),('reader',[bin/'reader_check_time',peram,'5184','KEYS_ONLY','1'])]
status_path=Path(str(point)+'.pipeline.json')
record={'status':'running','source_time':1,'sink_time':1,'physical_rhs':20736,'point_prefix':str(point),'standard_perambulator':str(peram),'started_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'commands':{name:[str(x) for x in cmd] for name,cmd in commands}}
record['MetalClover_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip()
record['Chroma_commit']=subprocess.check_output(['git','rev-parse','HEAD'],cwd='/Users/admin/Downloads/chroma.metalclover',text=True).strip()
record['pipeline_script_SHA256']=hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
def save():
 tmp=Path(str(status_path)+'.tmp');tmp.write_text(json.dumps(record,indent=2)+'\n');tmp.replace(status_path)
child=None
requested=[]
def stop(sig,_):
 requested.append(sig)
 if child is not None and child.poll() is None:child.send_signal(signal.SIGTERM)
signal.signal(signal.SIGINT,stop);signal.signal(signal.SIGTERM,stop)
start=time.monotonic()
for stage,cmd in commands:
 if requested:
  record.update(status='paused',stage=stage);save();raise SystemExit(75)
 record.update(stage=stage);save()
 with (scratch/('full-'+stage+'-controller.log')).open('x') as log:
  child=subprocess.Popen([str(x) for x in cmd],stdout=log,stderr=subprocess.STDOUT)
  record['controller_PID']=child.pid;save();code=child.wait()
 if code:
  record.update(status='paused' if code==75 else 'failed',exit_code=code,elapsed_seconds=time.monotonic()-start);save();raise SystemExit(code)
 if stage=='solve':
  completed=json.loads(Path(str(point)+'.complete.json').read_text())
  assert completed['physical_rhs']==20736 and completed['source_time']==completed['sink_time']==1
  assert completed['MG_hierarchy_builds_across_sessions']==completed['session_count']==1
 if stage=='project':
  completed=json.loads(Path(str(peram)+'.complete.json').read_text())
  assert completed['source_time']==completed['sink_time']==1 and completed['spin_records']==16
record.update(status='complete',stage='complete',elapsed_seconds=time.monotonic()-start,finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat());save()
print(json.dumps(record,indent=2),flush=True)
