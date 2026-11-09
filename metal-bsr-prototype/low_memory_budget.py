#!/usr/bin/env python3
"""Array accounting for LDL clover plus streamed CPU FP64 refinement.

Sparse graphs use the 33-neighbor Schur stencil bound. The measured Metal
allocation is still authoritative; this accounts for persistent buffers,
not driver/compiler/OS memory or transient CPU hierarchy construction.
"""
import argparse
import json
from math import prod,ceil
from pathlib import Path
p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--shape',type=int,nargs=4,required=True)
p.add_argument('--blocks',type=int,nargs=4,required=True)
p.add_argument('--restart',type=int,default=8)
p.add_argument('--out',type=Path,required=True)
a=p.parse_args();v=prod(a.shape);half=v//2;n=half*12;nv=8;bs=16
assert all(s%b==0 and (s//b)%2==0 for s,b in zip(a.shape,a.blocks))
rows1=prod(s//b for s,b in zip(a.shape,a.blocks));rows2=rows1//16;n1=rows1*bs;n2=rows2*bs
f={}
f['fine_weighted_gauge']=v*36*8
f['fine_even_and_odd_LDL_shared_with_smoother']=v*72*4
f['fine_neighbors_and_parity_maps']=v*36
f['fine_Schur_scratch']=2*n*8
f['fine_transfer_basis']=n*nv*8
f['fine_cycle_scratch']=3*n*8
f['fine_cycle_maps']=(2*half+rows1+1)*4
f['fine_cycle_coarse_rhs_and_solution']=2*n1*8
f['fine_FGMRES_vectors']=(2*a.restart+7)*n*8
f['fine_FGMRES_reductions']=a.restart*ceil(n/256)*8+ceil(n/256)*4+1024
f['reduction_and_reconstruction']=2*v*12*8

def matrix(label,rows):
 c=rows*bs;parts=ceil(c/256)
 f[label+'_matrix']=rows*33*bs*bs*8
 f[label+'_diagonal']=rows*bs*bs*8
 f[label+'_csr']=(rows+1+rows*33)*4
 f[label+'_eight_vectors']=8*c*8
 f[label+'_reductions']=2*parts*8+88
matrix('coarse1',rows1);matrix('bottom',rows2)
f['coarse1_FGMRES_vectors']=23*n1*8
f['coarse1_FGMRES_reductions']=8*ceil(n1/256)*8+ceil(n1/256)*4+1024
f['coarse1_cycle_transfer']=rows1*bs*nv*8
f['coarse1_cycle_maps']=(2*rows1+rows2+1)*4
f['coarse1_cycle_diagonal']=rows1*bs*bs*8
f['coarse1_cycle_scratch']=3*n1*8
f['coarse1_cycle_bottom_rhs_and_solution']=2*n2*8
slab=prod(a.shape[:3])*8
r={'shape':a.shape,'blocks':a.blocks,'fine_restart':a.restart,'vectors_per_chirality':nv,'coarse_dofs':[n1,n2],
 'persistent_GPU_bytes_by_group':f,'persistent_GPU_sum_GiB':sum(f.values())/2**30,
 'GPU_sum_plus_128MiB_padding_GiB':sum(f.values())/2**30+.125,
 'CPU_full_FP64_b_residual_clover_arrays':0,'CPU_FP64_x_file_GiB':v*12*16/2**30,
 'CPU_FP64_residual_slab_MiB':slab*12*16/2**20,
 'note':__doc__}
a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps({k:x for k,x in r.items() if k not in ['persistent_GPU_bytes_by_group','note']},indent=2))
