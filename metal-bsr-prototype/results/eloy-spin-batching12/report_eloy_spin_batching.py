from pathlib import Path
import csv,hashlib,json,math,re,subprocess
root=Path('/private/tmp/chroma-eloy-spin-batching')
chroma=Path('/Users/admin/Downloads/chroma.metalclover')
status=json.loads((root/'validation-status.json').read_text())
assert status['status']=='complete' and all(x['exit']==0 for x in status['steps'])
reg=(root/'regression.log').read_text()
assert 'PASS MetalBiCGClover and MetalMGClover:' in reg
fields=lambda s:dict(re.findall(r'(\w+)=([^ ]+)',s))
spin=[fields(line) for line in reg.splitlines() if line.startswith('ELOY_SPIN ')]
expected={(m,n,v,s) for m in ['MetalBiCGClover','MetalMGClover'] for n in [1,3] for v in range(n) for s in range(4)}
actual={(x['solver'],int(x['num_vecs']),int(x['vector']),int(x['source_spin'])) for x in spin}
assert len(spin)==32 and actual==expected
for x in spin:
 assert int(x['physical_rhs'])==4*int(x['num_vecs'])
 assert int(x['max_rhs'])==(1 if int(x['num_vecs'])==1 else 2)
 assert 0<=float(x['full_relative'])<=1e-12
 assert 0<=float(x['projected_error'])<=1e-11
batches=[fields(line) for line in reg.splitlines() if line.startswith('ELOY_BATCH ')]
assert len(batches)==4
for x in batches:
 n=int(x['num_vecs']);widths=[int(k) for k in x['widths'].split(',')]
 assert widths==([1,1,1,1] if n==1 else [2,1,2,1,2,1,2,1])
 assert int(x['qprop_calls'])==len(widths) and sum(widths)==4*n
native=[line for line in reg.splitlines() if line.startswith(('HARNESS ','RECONSTRUCTED ','QPROP_ARRAY ','SUPERB_ARRAY ','ELOY_SPIN '))]
assert len(native)==64
smoke=(root/'spin-major-smoke.log').read_text()
assert 'PASS real 12^3x96 Eloy spin-major distillation qprop' in smoke
assert 'submission_order=spin_major' in smoke and 't_source=0' in smoke
assert smoke.count('NEAR_NULL_SVD ')==1
with (root/'spins1.csv').open() as f:rows=list(csv.DictReader(f))
expected={(m,n,c) for m in ['MetalBiCGClover','MetalMGClover'] for n,c in [(1,1),(1,8),(1,16),(2,16)]}
assert len(rows)==8 and {(x['method'],int(x['num_vecs']),int(x['caller_tile'])) for x in rows}==expected
observed={};current=None
for line in smoke.splitlines():
 x=fields(line)
 if line.startswith('BATCH_BEGIN '):
  current=(x['method'],int(x['num_vecs']),int(x['caller_tile']))
  observed[current]={'qprop_calls':0,'GPU_tiles':0,'maximum_GPU_rhs':0}
 elif ' array: rhs=' in line:
  assert current is not None
  out=observed[current];out['qprop_calls']+=1;out['GPU_tiles']+=int(x['gpu_tiles'])
  out['maximum_GPU_rhs']=max(out['maximum_GPU_rhs'],int(x['resident_GPU_rhs']))
 elif line.startswith('BATCH_RESULT '):current=None
for x in rows:
 n=int(x['num_vecs']);cap=int(x['caller_tile']);rhs=int(x['rhs'])
 assert int(x['source_spins'])==4 and rhs==4*n
 assert 0<=float(x['max_full_relative'])<=1e-12 and 0<=float(x['max_solution_difference'])<=1e-8
 assert math.isclose(float(x['seconds_per_colorvec']),4*float(x['seconds_per_rhs']),rel_tol=1e-14)
 out=observed[(x['method'],n,cap)]
 assert out['qprop_calls']==out['GPU_tiles']==4*math.ceil(n/cap)
 assert out['maximum_GPU_rhs']==min(n,cap)
 x.update(out)
manifest=json.loads((root/'sources.json').read_text())
for name,record in manifest.items():
 p=chroma/name
 assert {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}==record
 assert p.read_bytes()==(Path('/private/tmp/chroma-metal-integration')/name).read_bytes()
restored=['lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.cc','lib/meas/inline/hadron/inline_prop_and_matelem_distillation_superb_w.h','lib/util/ferm/mgproton.h']
for name in restored:
 assert (chroma/name).read_bytes()==subprocess.run(['git','show','b3a380ea6:'+name],cwd=chroma,check=True,capture_output=True).stdout
summary={'status':'pass','lattice':[12,12,12,96],'mass':-0.0743,'source_time':0,'source_spins':4,'submission_order':'source spin outer, eigenvector batches inner','max_rhs_rounding':False,'max_rhs_1_verified':True,'native_checks':len(native),'spin_source_checks':len(spin),'max_full_relative':max(float(x['full_relative']) for x in spin),'max_projected_element_error':max(float(x['projected_error']) for x in spin),'smoke_max_full_relative':max(float(x['max_full_relative']) for x in rows),'smoke_max_solution_difference':max(float(x['max_solution_difference']) for x in rows),'one_reused_MG_hierarchy_per_smoke_solver':True,'byte_exact_restoration_from':'b3a380ea62c7c9002781a5bc1e9a557e4d3f3521','byte_exact_restored_files':restored,'caller_batches_verified':batches,'spin_source_checks_detail':spin,'smoke_measurements':rows,'scope':'Real-gauge qprop/colorvec source assembly and direct contraction checks; full inline perambulator file I/O and full 5184-vector run were not run.'}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['caller_batches_verified','spin_source_checks_detail','smoke_measurements']},indent=2))
