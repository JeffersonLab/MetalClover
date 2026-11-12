#!/usr/bin/env python3
"""Embed the runtime Metal library so an installed dylib needs no source tree."""
from pathlib import Path
import sys
names = ('kernels.metal','bicgstab.metal','clover.metal','multigrid.metal',
         'sparse_coarse.metal','compact_clover.metal','mg_real_tune.metal',
         'packed_clover.metal')
root = Path(__file__).resolve().parent
source = '\n'.join((root / name).read_text() for name in names)
assert ')MC_SHADER"' not in source
Path(sys.argv[1]).write_text('namespace MetalCloverEnv {\nstatic const char metalclover_shader[] = R"MC_SHADER(\n'
                           + source + ')MC_SHADER";\n} // namespace MetalCloverEnv\n')
