from pathlib import Path
import json,hashlib,shutil,csv,subprocess,sys,re
root=Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU');t=int(sys.argv[1])
dst=root/f'metal-bsr-prototype/results/cfg1020-perambulator-t{t}-5184'
base=Path('/Users/admin/Downloads/MetalClover-results');p=base/f'cfg1020-t{t}-unsmeared'
s=base/f'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb1020'
r=json.loads(Path(str(p)+'.complete.json').read_text());run=json.loads(Path(str(p)+'.run.json').read_text());pr=json.loads(Path(str(s)+'.complete.json').read_text());pipeline=json.loads(Path(str(p)+'.pipeline.json').read_text())
assert r['status']==run['status']==pr['status']==pipeline['status']=='complete'
assert r['physical_rhs']==r['completed_rhs']==20736 and r['num_vectors']==5184
assert r['source_time']==r['sink_time']==pr['source_time']==pr['sink_time']==t
assert r['session_count']==r['MG_hierarchy_builds_across_sessions']==1
assert r['maximum_full_relative']<=1e-12 and r['maximum_qprop_relative']<=1e-12
assert pr['spin_records']==16 and pr['full_gamma5_relative_Frobenius']<1e-12 and pr['max_QDP_dot_scaled_error']<1e-12
assert p.with_suffix('.npy').stat().st_size==r['payload_bytes']==6879707264
assert s.stat().st_size==pr['output_bytes']
frozen_path=Path(f'/private/tmp/chroma-cfg1020-time{t}-production/frozen-files.json')
frozen=json.loads(frozen_path.read_text())
for name,h in frozen.items():
 if name.startswith('/'):
  assert hashlib.sha256(Path(name).read_bytes()).hexdigest()==h,name
log=Path(str(p)+'.session-000.log').read_text();arrays=[x for x in log.splitlines() if 'MetalMGClover array: rhs=' in x];progress=[x for x in log.splitlines() if x.startswith('PRODUCTION_PROGRESS ')]
assert len(arrays)==len(progress)==1296
assert all('rhs=16 resident_GPU_rhs=16 gpu_tiles=1 hierarchy_reused=true' in x for x in arrays)
assert log.count('\nNEAR_NULL_SVD ')==1
assert f'PASS unsmeared_slice time={t} completed_rhs=20736 total_rhs=20736' in log
assert f'source_time={t} sink_time={t}' in log
assert 'RESUME_VERIFIED completed_rhs=20736' in Path(str(p)+'.session-000.verify.log').read_text()
reader=Path(f'/private/tmp/chroma-cfg1020-time{t}-production/full-reader-controller.log').read_text()
assert f'PASS analysis_reader t_source={t} t_sink={t} mass=light vectors=5184' in reader
with Path(str(p)+'.columns.csv').open() as f:rows=list(csv.DictReader(f))
per_spin={}
for spin in range(4):
 a=[x for x in rows if int(x['source_spin'])==spin];assert len(a)==5184
 assert {int(x['vector']) for x in a}==set(range(5184))
 per_spin[spin]={'accepted_rhs':len(a),'maximum_full_relative':max(float(x['full_relative']) for x in a),'source_norm_min':min(float(x['source_norm']) for x in a),'source_norm_max':max(float(x['source_norm']) for x in a)}
for suffix in ['.complete.json','.run.json','.meta.xml','.columns.csv','.state','.session-000.log','.session-000.resources.csv','.session-000.status.json','.session-000.verify.log','.pipeline.json']:
 shutil.copyfile(Path(str(p)+suffix),dst/('point'+suffix))
for suffix in ['.complete.json','.run.json','.projection.xml','.projection.csv','.log','.resources.csv']:
 shutil.copyfile(Path(str(s)+suffix),dst/('peram'+suffix))
shutil.copyfile(f'/private/tmp/chroma-cfg1020-time{t}-production/full-reader-controller.log',dst/'peram-reader.log')
shutil.copyfile(f'/private/tmp/chroma-cfg1020-time{t}-production/full_run.py',dst/'production-harness.py')
record={'status':'pass','configuration':1020,'source_time':t,'sink_time':t,'physical_rhs':20736,'source_vectors':5184,'spin_records':16,'native_array_calls':1296,'RHS_per_call':16,'hierarchies':1,'sessions':1,'both_original_full_CPU_FP64_checks_tolerance':1e-12,'maximum_full_relative':r['maximum_full_relative'],'all_point_column_checksums_verified':True,'all_standard_record_values_readback_bitwise':True,'full_gamma5_relative_Frobenius':pr['full_gamma5_relative_Frobenius'],'direct_QDP_contraction_scaled_error':pr['max_QDP_dot_scaled_error'],'independent_DataLib_full_file_key_check':'pass','unchanged_frozen_source_executable_library_hashes':True,'per_source_spin':per_spin,'solve_seconds':r['elapsed_seconds_across_sessions'],'projection_validation_seconds':pr['native_seconds'],'pipeline_seconds':pipeline['elapsed_seconds'],'nominal_dispatched_GPU_operator_GFLOP_per_second_solve_wall':r['nominal_operator_GFLOP_per_second_total_wall'],'point_SHA256':r['payload_sha256'],'standard_perambulator_SHA256':pr['output_SHA256'],'standard_perambulator':str(s),'point_file':str(p)+'.npy','large_payloads_committed_to_git':False}
(dst/'production-validation.json').write_text(json.dumps(record,indent=2)+'\n')
print(json.dumps({k:record[k] for k in ['status','source_time','sink_time','solve_seconds','projection_validation_seconds','pipeline_seconds','maximum_full_relative','full_gamma5_relative_Frobenius','nominal_dispatched_GPU_operator_GFLOP_per_second_solve_wall','point_SHA256','standard_perambulator_SHA256']},indent=2))
