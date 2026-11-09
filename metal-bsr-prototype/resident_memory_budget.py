#!/usr/bin/env python3
"""Conservative lower bound for this resident three-level Metal solver.

Only primary live arrays are counted: excludes sparse coarse matrices,
coarse solver workspaces, host geometry maps, Metal heaps/driver overhead,
allocator padding, transient hierarchy construction, and the OS/apps.
"""
import argparse
import json
from math import prod
from pathlib import Path

p = argparse.ArgumentParser(description=__doc__)
p.add_argument('--shape', type=int, nargs=4, required=True)
p.add_argument('--restart', type=int, default=8)
p.add_argument('--vectors-per-chirality', type=int, default=12)
p.add_argument('--resident-rhs', type=int, default=1)
p.add_argument('--out', type=Path, required=True)
p.add_argument('--device-limits', type=Path, help='JSON from metal_memory_limits')
a = p.parse_args()
assert min(a.shape)>0 and all(n%2==0 for n in a.shape)
assert a.restart>0 and a.vectors_per_chirality>0 and a.resident_rhs>0
v = prod(a.shape); half = v//2; r = a.resident_rhs
fields = {
    'host_fp64_physical_gauge': v*4*9*16,
    'host_fp64_even_and_odd_clover': v*2*6*6*16,
    'host_fp64_b_x_residual': v*12*16*3*r,
    'metal_fp32_weighted_gauge': v*4*9*8,
    'metal_fp32_even_clover_inverse': half*2*6*6*8,
    'metal_fp32_odd_clover': half*2*6*6*8,
    'metal_neighbors_and_parity_maps': v*8*4+v*4,
    'metal_schur_two_half_fields': half*12*8*2*r,
    'metal_fine_smoother_inverse': half*2*6*6*8,
    'metal_fine_transfer_basis': half*12*a.vectors_per_chirality*8,
    'metal_fine_cycle_three_half_fields': half*12*8*3*r,
    'metal_fine_fgmres_fields': half*12*8*(6+(a.restart+1)+a.restart)*r,
    'metal_correction_full_and_two_half_fields': (v+2*half)*12*8*r,
}
total = sum(fields.values())
known_buffers = {
    'fine_transfer_basis': half*12*a.vectors_per_chirality*8,
    'fine_fgmres_V': half*12*8*(a.restart+1)*r,
    'fine_fgmres_Z': half*12*8*a.restart*r,
    'weighted_gauge': v*4*9*8,
    'even_clover_inverse': half*2*6*6*8,
    'correction_full_field': v*12*8*r,
    'neighbors': v*8*4,
}
result = {
    'shape': a.shape,
    'volume': v,
    'resident_rhs': r,
    'fine_restart': a.restart,
    'vectors_per_chirality': a.vectors_per_chirality,
    'bytes_by_live_array_group': fields,
    'gib_by_live_array_group': {k:n/2**30 for k,n in fields.items()},
    'lower_bound_bytes': total,
    'lower_bound_gib': total/2**30,
    'known_single_metal_buffers_bytes': known_buffers,
    'largest_known_single_metal_buffer_bytes': max(known_buffers.values()),
    'excludes': __doc__.split('Only primary live arrays are counted: ')[1].strip(),
    'fits_32_gib_even_before_excluded_costs': total<=32*2**30,
    'source': 'compact_clover.hpp, compact_solve.hpp, fgmres.hpp, three_level.hpp',
}
if a.device_limits:
    device = json.loads(a.device_limits.read_text())
    result['device_limits'] = device
    result['known_buffers_within_device_max_buffer_length'] = max(known_buffers.values())<=device['max_buffer_length_bytes']
    result['live_storage_within_physical_memory_before_excluded_costs'] = total<=device['physical_memory_bytes']
a.out.parent.mkdir(parents=True, exist_ok=True)
a.out.write_text(json.dumps(result, indent=2)+'\n')
print(json.dumps(result, indent=2))
