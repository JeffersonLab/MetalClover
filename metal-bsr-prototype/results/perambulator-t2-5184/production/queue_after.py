from pathlib import Path
import datetime,hashlib,json,os,signal,subprocess,time
root=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
work=Path('/private/tmp/chroma-time2-production')
point1=Path('/Users/admin/Downloads/MetalClover-results/cfg1000-t1-unsmeared')
point2=Path('/Users/admin/Downloads/MetalClover-results/cfg1000-t2-unsmeared')
status=Path(str(point2)+'.queue.json')
if status.exists():raise RuntimeError('time-two queue already exists')
def read(path):return json.loads(Path(path).read_text())
def sha(path):
 h=hashlib.sha256()
 with Path(path).open('rb') as f:
  for b in iter(lambda:f.read(4*1024**2),b''):h.update(b)
 return h.hexdigest()
record={'status':'queued','source_time':2,'sink_time':2,'physical_rhs':20736,'queue_PID':os.getpid(),'dependency_pipeline':str(point1)+'.pipeline.json','pipeline_script':str(work/'full_run.py'),'queued_UTC':datetime.datetime.now(datetime.timezone.utc).isoformat(),'continues_without_iPhone_connection':True}
def save():
 temp=Path(str(status)+'.tmp');temp.write_text(json.dumps(record,indent=2)+'\n');temp.replace(status)
stop=[];child=None

def on_signal(sig,_):
 stop.append(sig)
 if child is not None and child.poll() is None:child.send_signal(signal.SIGTERM)
signal.signal(signal.SIGTERM,on_signal);signal.signal(signal.SIGINT,on_signal)
keeper=subprocess.Popen(['/usr/bin/caffeinate','-i','-w',str(os.getpid())],stdin=subprocess.DEVNULL,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
record['awake_keeper_PID']=keeper.pid;save()
try:
 while not stop:
  dependency=read(str(point1)+'.pipeline.json')
  record['dependency_status']=dependency.get('status');record['dependency_stage']=dependency.get('stage')
  record['last_checked_UTC']=datetime.datetime.now(datetime.timezone.utc).isoformat();save()
  if dependency['status']=='complete':break
  time.sleep(30)
 if stop:
  record['status']='paused';save();raise SystemExit(75)
 # Completion includes the standard-file verification and independent reader.
 first=read(str(point1)+'.complete.json')
 peram1=Path(dependency['standard_perambulator']);projected=read(str(peram1)+'.complete.json')
 assert first['source_time']==first['sink_time']==1 and first['completed_rhs']==20736
 assert first['maximum_full_relative']<=1e-12 and first['maximum_qprop_relative']<=1e-12
 assert projected['source_time']==projected['sink_time']==1 and projected['spin_records']==16
 assert projected['status']=='complete' and peram1.stat().st_size==projected['output_bytes']
 staged=read(root/'metal-bsr-prototype/results/perambulator-t2-5184/staged-run.json')
 assert sha(staged['eigenvector_file'])==staged['eigenvector_SHA256']
 assert sha(staged['gauge_file'])==staged['gauge_SHA256']
 for name,value in read(root/'metal-bsr-prototype/results/perambulator-t2-5184/frozen-files.json').items():
  if name.startswith('/'):assert sha(name)==value,name
 record.update(status='running',dependency_standard_perambulator_SHA256=projected['output_SHA256'],started_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat())
 with (work/'pipeline-controller.log').open('x') as log:
  child=subprocess.Popen(['python3',str(work/'full_run.py')],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT)
  record['pipeline_PID']=child.pid;save();code=child.wait()
 record.update(status='complete' if code==0 else ('paused' if code==75 else 'failed'),exit_code=code,finished_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat());save()
 raise SystemExit(code)
except Exception as error:
 record.update(status='failed',error=str(error));save();raise
finally:
 if keeper.poll() is None:keeper.terminate()
 keeper.wait()
