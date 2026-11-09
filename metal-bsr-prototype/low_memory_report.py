#!/usr/bin/env python3
"""Summarize only completed, original-Chroma-verified 32³ GPU solves."""
import argparse
import csv
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import re
p=argparse.ArgumentParser(description=__doc__);p.add_argument('directory',type=Path);a=p.parse_args();d=a.directory

def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))

settings=json.loads((d/'settings.json').read_text());input_data=json.loads((d/'input.json').read_text())
result={}
for short,method in [('bicg','metal-bicgstab'),('mg','metal-mg')]:
    measured=rows(d/f'test.{short}.csv');native=rows(d/f'test.{short}-native.csv');log=(d/f'{short}.log').read_text();native_log=(d/f'verify-{short}.log').read_text()
    if 'PASS all full FP64 residuals' not in log or 'PASS original Chroma FP64 verification' not in native_log:raise RuntimeError('incomplete '+short)
    for table in [measured,native]:
        if len(table)!=2 or sorted(int(r['rhs']) for r in table)!=[0,1]:raise RuntimeError('missing or duplicate fields')
        for r in table:
            rho=float(r['full_relative'])
            if r['method']!=method or not math.isfinite(rho) or not 0<=rho<=1e-12:raise RuntimeError('full FP64 tolerance rejected')
    hierarchy=0;filter_calls=0
    if short=='mg':
        if log.count('MG_READY ')!=1:raise RuntimeError('expected exactly one reused hierarchy')
        ready=re.search(r'MG_READY hierarchy_seconds=([\d.eE+-]+) hierarchy_fine_S=(\d+) GPU_GiB=([\d.eE+-]+) recommended_GiB=([\d.eE+-]+)',log)
        hierarchy=float(ready[1]);filter_calls=int(ready[2]);gpu=float(ready[3]);recommended=float(ready[4])
    else:
        gpu=float(re.search(r'BICG_READY GPU_GiB=([\d.eE+-]+)',log)[1]);recommended=24.960006713867188
    elapsed=sum(float(r['seconds']) for r in measured)
    result[method]={'rhs_count':2,'mean_seconds_per_rhs':elapsed/2,'solve_seconds':elapsed,'hierarchy_seconds':hierarchy,'mean_seconds_per_rhs_with_hierarchy':(elapsed+hierarchy)/2,'hierarchy_fine_S':filter_calls,'fine_S_total':sum(int(r['fine_schur_applications']) for r in measured),'krylov_total':sum(int(r['krylov']) for r in measured),'outer_corrections':[int(r['outer']) for r in measured],'CPU_FP64_residual_seconds':sum(float(r['residual_seconds']) for r in measured),'max_original_Chroma_FP64_relative':max(float(r['full_relative']) for r in native),'GPU_allocated_GiB_at_ready':gpu,'GPU_recommended_working_set_GiB':recommended,'standalone_original_FP64_probe':float(re.search(r'FP64 original_Chroma_operator_match=([\d.eE+-]+)',log)[1]),'measured':measured,'native':native}
bi=result['metal-bicgstab'];mg=result['metal-mg'];input_data['execution_status']='Completed; both solvers and all four original Chroma FP64 checks passed at1e-12'
settings['actual_mg']={'fine_aggregation':[4,4,4,4],'vectors_per_chirality_each_level':8,'fine_GPU_filter_sweeps':256,'second_FP64_filter_sweeps':64,'fine_restart':8,'pre_steps':1,'post_steps':1,'smoother_weight':1,'coarse_restart':8,'coarse_limit':32,'bottom_limit':8,'coarse_and_bottom_tolerance':.1,'hierarchies_built':1,'resident_rhs':1}
summary={'completed_utc':datetime.now(timezone.utc).isoformat(),'input':input_data,'settings':settings,'summary':result,'comparisons':{'BiCGStab_over_MG_solve_time':bi['solve_seconds']/mg['solve_seconds'],'BiCGStab_over_MG_with_hierarchy':bi['mean_seconds_per_rhs']/mg['mean_seconds_per_rhs_with_hierarchy'],'fine_S_reduction_vs_BiCGStab':1-mg['fine_S_total']/bi['fine_S_total']},'precision':'GPU FP32 LDL clover and Krylov work; original-gauge CPU FP64 clover recomputation, accumulation, full residual and norm; final original Chroma UnprecCloverLinOp FP64 verification.','source_tags':[20261107,20261108],'source_initial_guess':'zero for every solver/RHS','data_files':json.loads((d/'accepted-files.json').read_text())}
(d/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');(d/'input.json').write_text(json.dumps(input_data,indent=2)+'\n')
print(json.dumps({'summary':result,'comparisons':summary['comparisons']},indent=2))
