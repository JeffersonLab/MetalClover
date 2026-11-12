from pathlib import Path
import os,json,subprocess,time
root=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU/metal-bsr-prototype')
out=Path('/private/tmp/mg-null32');out.mkdir(exist_ok=True)
xml=subprocess.check_output(['git','show','HEAD:metal-bsr-prototype/results/performance32-fit/action.xml'],cwd=root.parent)
(out/'action.xml').write_bytes(xml)
state={'status':'running','started':time.time(),'source_tags':[20261107,20261108],'steps':[]}
def save():
 p=out/'state.next';p.write_text(json.dumps(state,indent=2)+'\n');p.replace(out/'state.json')
def run(name,args):
 state['phase']=name;save();begin=time.time()
 code=subprocess.call(args,cwd='/private/tmp/mg-null32-runtime',env={**os.environ,'VECLIB_MAXIMUM_THREADS':'1'})
 state['steps'].append({'name':name,'exit':code,'elapsed':time.time()-begin});save()
 if code:raise RuntimeError(f'{name}: exit {code}')
guard=['python3',str(root/'run_guarded.py')]
try:
 run('mg',guard+['--log',str(out/'mg.log'),'--resources',str(out/'mg.resources.csv'),'--','/private/tmp/mg-null32-solve','/private/tmp/mg-fit32',str(out/'near'),'mg','2','4,4,4,4','8','256','64','8','--basis','near-null-svd','--null-vectors','24','--null-iters','128','--null-restart','8','--null-warm','32','--null-tol','.001','--null-batch','4'])
 run('mg-native',guard+['--log',str(out/'mg.native.log'),'--resources',str(out/'mg.native.resources.csv'),'--','/private/tmp/mg-near-null-verify',str(out/'action.xml'),'/private/tmp/mg-fit32',str(out/'near'),'metal-mg','2',str(out/'mg.native.csv')])
 run('bicg',guard+['--log',str(out/'bicg.log'),'--resources',str(out/'bicg.resources.csv'),'--','/private/tmp/mg-null32-solve','/private/tmp/mg-fit32',str(out/'bicg'),'bicg','2','4,4,4,4','8','256','64','8'])
 run('bicg-native',guard+['--log',str(out/'bicg.native.log'),'--resources',str(out/'bicg.native.resources.csv'),'--','/private/tmp/mg-near-null-verify',str(out/'action.xml'),'/private/tmp/mg-fit32',str(out/'bicg'),'metal-bicgstab','2',str(out/'bicg.native.csv')])
 state['status']='complete'
except BaseException as e:
 state['status']='failed';state['error']=str(e);raise
finally:state['finished']=time.time();save()
