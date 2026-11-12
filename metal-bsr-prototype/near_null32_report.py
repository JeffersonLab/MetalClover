#!/usr/bin/env python3
"""Validate two-source 32-cubed comparison with scoped operator-work rates."""
import csv,json,math,re,statistics,sys
from pathlib import Path
root=Path(sys.argv[1])
def rows(name):
 with (root/name).open() as f:return list(csv.DictReader(f))
def number(log,key):
 value=float(re.search(r'\b'+key+r'=([0-9.eE+-]+)',log)[1])
 if not math.isfinite(value):raise ValueError(key)
 return value
methods={}
for name,table,logname,nativefile,method in (
 ('near-null-svd','near.mg.csv','mg.log','mg.native.csv','metal-mg'),
 ('metal-bicgstab','bicg.bicg.csv','bicg.log','bicg.native.csv','metal-bicgstab')):
 data=rows(table);native=rows(nativefile);log=(root/logname).read_text()
 if len(data)!=2 or len(native)!=2 or {int(x['rhs']) for x in data}!={0,1} or {int(x['rhs']) for x in native}!={0,1}:
  raise ValueError('missing/duplicated RHS')
 for x in data+native:
  residual=float(x['full_relative'])
  if x['method']!=method or not math.isfinite(residual) or not 0<=residual<=1e-12:raise ValueError('FP64 acceptance')
 if 'PASS all full FP64 residuals' not in log:raise ValueError('unfinished solve')
 if number(log,'original_Chroma_operator_match')>1e-13:raise ValueError('operator mismatch')
 for x in data:
  flops=int(x['nominal_operator_flops'])
  parts=sum(int(x[k]) for k in ('fine_flops','coarse_flops','transfer_flops','local_vector_flops'))
  if flops!=parts or abs(float(x['operator_gflops'])*float(x['seconds'])*1e9/flops-1)>1e-12:raise ValueError('work accounting')
 setup=number(log,'hierarchy_seconds') if method=='metal-mg' else 0
 if method=='metal-mg' and log.count('MG_READY ')!=1:raise ValueError('hierarchy rebuilt')
 mean=statistics.fmean(float(x['seconds']) for x in data)
 methods[name]={
  'seconds':[float(x['seconds']) for x in data], 'mean_seconds_per_rhs':mean,
  'hierarchy_seconds':setup,'measured_seconds_per_rhs_including_setup':mean+setup/2,
  'krylov_iterations':[int(x['krylov']) for x in data],
  'fine_schur_applications':[int(x['fine_schur_applications']) for x in data],
  'CPU_FP64_residual_seconds':[float(x['residual_seconds']) for x in data],
  'maximum_original_Chroma_residual':max(float(x['full_relative']) for x in native),
  'maximum_inline_residual':max(float(x['full_relative']) for x in data),
  'nominal_operator_flops':[int(x['nominal_operator_flops']) for x in data],
  'nominal_operator_GFLOP_per_second':[float(x['operator_gflops']) for x in data],
  'weighted_operator_GFLOP_per_second':sum(int(x['nominal_operator_flops']) for x in data)/(2*mean*1e9),
  'GPU_allocation_at_ready_GiB':number(log,'GPU_GiB'),
  'peak_RSS_GiB':number(log,'peak_RSS_GiB'),
  'projected_seconds_per_rhs_including_setup':{str(n):mean+setup/n for n in (96,192)}}
near=methods['near-null-svd'];bi=methods['metal-bicgstab']
summary={'status':'complete','shape':[32,32,32,256],'mass':-.0860,'measured_rhs':2,
 'rhs_tags':[20261107,20261108],'resident_rhs':1,'hierarchies':1,'methods':methods,
 'warm_speedup':bi['mean_seconds_per_rhs']/near['mean_seconds_per_rhs'],
 'speedup_including_measured_setup':bi['mean_seconds_per_rhs']/near['measured_seconds_per_rhs_including_setup'],
 'fine_work_reduction':1-sum(near['fine_schur_applications'])/sum(bi['fine_schur_applications']),
 'rate_scope':'nominal dispatched operator work / end-to-end solve seconds; see OPERATOR_WORK.md',
 'basis_svd':{k:float(v) for k,v in re.findall(r'(\w+)=([\d.eE+-]+)',next(l for l in (root/'mg.log').read_text().splitlines() if l.startswith('NEAR_NULL_SVD ')))}}
(root/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
