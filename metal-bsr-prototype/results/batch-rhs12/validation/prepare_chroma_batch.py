from pathlib import Path
import subprocess,shutil,os,json,shlex
source=Path('/Users/admin/Downloads/chroma.metalclover')
old=Path('/private/tmp/chroma-metal-integration')
out=Path('/private/tmp/chroma-batch-validation');out.mkdir(exist_ok=True)
# Preserve the old installed binary and use the already built matching Chroma
# archive as an incremental build. The measured code is copied byte for byte.
oldprefix=Path('/private/tmp/metalclover-install')
backup=out/'previous-libmetalclover.dylib'
if not backup.exists():shutil.copyfile(oldprefix/'lib/libmetalclover.dylib',backup)
shutil.copyfile('/private/tmp/metalclover-batch-install/lib/libmetalclover.dylib',oldprefix/'lib/libmetalclover.dylib')
shutil.copyfile('/private/tmp/metalclover-batch-install/include/metalclover/metalclover.h',oldprefix/'include/metalclover/metalclover.h')
files=['lib/actions/ferm/invert/metal_solvers/syssolver_linop_metal_clover.h','lib/actions/ferm/invert/metal_solvers/syssolver_linop_metal_clover.cc','mainprogs/tests/t_metal_clover.cc','mainprogs/tests/t_metal_clover_batch.cc']
for name in files:shutil.copyfile(source/name,old/name)
b=Path('/private/tmp/chroma-metal-validation/enabled/mainprogs/tests')
(b/'Makefile.batch').write_text('''include Makefile

t_metal_clover_batch.o: /private/tmp/chroma-metal-integration/mainprogs/tests/t_metal_clover_batch.cc
\t$(CXX) $(DEFS) $(DEFAULT_INCLUDES) $(INCLUDES) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) -c -o $@ $<
t_metal_clover_batch: t_metal_clover_batch.o $(t_metal_clover_DEPENDENCIES)
\t$(CXXLINK) t_metal_clover_batch.o $(LDADD) $(LIBS)
''')
for name,command in [('regression-build',['make','t_metal_clover']),('benchmark-build',['make','-f','Makefile.batch','t_metal_clover_batch'])]:
 print(name,flush=True)
 with (out/(name+'.log')).open('w') as log:
  subprocess.run(command,cwd=b,stdout=log,stderr=subprocess.STDOUT,check=True)
print('NATIVE_BATCH_BUILD_COMPLETE',flush=True)
