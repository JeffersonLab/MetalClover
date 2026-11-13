from pathlib import Path
import csv,hashlib,json,math,re,shutil
root=Path('/private/tmp/chroma-distillation-spin')
repo=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
chroma=Path('/Users/admin/Downloads/chroma.metalclover')
status=json.loads((root/'validation-status.json').read_text())
assert status['status']=='complete' and all(s['exit']==0 for s in status['steps'])
reg=(root/'regression.log').read_text()
assert 'PASS MetalBiCGClover and MetalMGClover:' in reg
fields=lambda line:dict(re.findall(r'(\w+)=([^ ]+)',line))
spin=[fields(line) for line in reg.splitlines() if line.startswith('DISTILLATION_SPIN ')]
expected={(m,n,v,s) for m in ['MetalBiCGClover','MetalMGClover'] for n in [1,3] for v in range(n) for s in range(4)}
actual={(f['solver'],int(f['num_vecs']),int(f['vector']),int(f['source_spin'])) for f in spin}
assert len(spin)==32 and actual==expected
for f in spin:
 assert int(f['physical_rhs'])==4*int(f['num_vecs'])
 assert int(f['effective_max_rhs'])==(4 if int(f['num_vecs'])==1 else 8)
 assert 0<=float(f['full_relative'])<=1e-12
 assert 0<=float(f['projected_error'])<=1e-11
native=[]
for line in reg.splitlines():
 if line.startswith(('HARNESS ','RECONSTRUCTED ','QPROP_ARRAY ','SUPERB_ARRAY ','DISTILLATION_SPIN ')):
  f=fields(line)
  native.append({'kind':line.split()[0],**f})
assert len(native)==64
smoke=(root/'four-spin-smoke.log').read_text()
assert 'PASS real 12^3x96 four-spin distillation qprop' in smoke
assert smoke.count('NEAR_NULL_SVD ')==1
with (root/'spins1.csv').open() as f: rows=list(csv.DictReader(f))
expected={(m,n,t) for m in ['MetalBiCGClover','MetalMGClover'] for n,t in [(1,1),(1,8),(1,16),(2,16)]}
actual={(x['method'],int(x['num_vecs']),int(x['caller_tile'])) for x in rows}
assert len(rows)==8 and actual==expected
observed={};current=None
for line in smoke.splitlines():
 f=fields(line)
 if line.startswith('BATCH_BEGIN '):
  current=(f['method'],int(f['num_vecs']),int(f['caller_tile']))
  observed[current]={'qprop_calls':0,'maximum_GPU_rhs':0,'GPU_tiles':0}
 elif ' array: rhs=' in line:
  assert current is not None
  item=observed[current];item['qprop_calls']+=1
  item['maximum_GPU_rhs']=max(item['maximum_GPU_rhs'],int(f['resident_GPU_rhs']))
  item['GPU_tiles']+=int(f['gpu_tiles'])
 elif line.startswith('BATCH_RESULT '): current=None
for row in rows:
 n=int(row['num_vecs']);rhs=int(row['rhs']);tile=int(row['caller_tile'])
 assert int(row['source_spins'])==4 and rhs==4*n
 assert 0<=float(row['max_full_relative'])<=1e-12
 assert 0<=float(row['max_solution_difference'])<=1e-8
 assert math.isclose(float(row['seconds_per_colorvec']),4*float(row['seconds_per_rhs']),rel_tol=1e-14)
 key=(row['method'],n,tile);item=observed[key]
 assert item['maximum_GPU_rhs']==min(rhs,tile)
 assert item['GPU_tiles']==item['qprop_calls']==math.ceil(rhs/tile)
 row.update(item)
manifest=json.loads((root/'sources.json').read_text())
for name,record in manifest.items():
 p=chroma/name
 assert {'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()}==record
 assert p.read_bytes()==(Path('/private/tmp/chroma-metal-integration')/name).read_bytes()
summary={'status':'pass','lattice':[12,12,12,96],'mass':-0.0743,'source_time':80,
 'source_spins':4,'rhs_index':'4 * vector + source_spin','tolerance':1e-12,
 'normalization':'positive cap rounded down to a multiple of four, at least four, bounded by total; nonpositive selects all',
 'native_checks':len(native),'four_spin_checks':len(spin),
 'four_spin_max_full_relative':max(float(f['full_relative']) for f in spin),
 'max_projected_element_error':max(float(f['projected_error']) for f in spin),
 'smoke_max_full_relative':max(float(f['max_full_relative']) for f in rows),
 'smoke_max_solution_difference':max(float(f['max_solution_difference']) for f in rows),
 'one_reused_MG_hierarchy_per_smoke_solver':True,
 'smoke_measurements':rows,'four_spin_native_checks':spin,
 'scope':'Real-gauge qprop/colorvec helper and tiled contractions; full inline file-I/O job not run. Small smoke timings are not 96/192-RHS performance measurements.'}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps({k:v for k,v in summary.items() if k not in ['smoke_measurements','four_spin_native_checks']},indent=2))
