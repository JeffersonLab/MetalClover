from pathlib import Path
import subprocess,shlex,json,os
root=Path('/Users/admin/Downloads/chroma.metalclover');out=Path('/private/tmp/chroma-batch-validation');out.mkdir(exist_ok=True)
with (out/'bootstrap.log').open('w') as log:
 subprocess.run(['autoreconf','-fi'],cwd=root,stdout=log,stderr=subprocess.STDOUT,check=True)
old=json.loads(Path('/private/tmp/chroma-metal-validation/validation.json').read_text())
args=shlex.split(old['configure_argv'])
args=[a for a in args if not a.startswith(('--with-metal=','--prefix='))]
for name,flag in [('enabled','--with-metal=/private/tmp/metalclover-batch-install'),('disabled','--without-metal')]:
 b=out/('configure-'+name);b.mkdir(exist_ok=True)
 print('configure-'+name,flush=True)
 with (out/('configure-'+name+'.log')).open('w') as log:
  subprocess.run([str(root/'configure'),*args,'--prefix='+str(out/'install'),flag],cwd=b,stdout=log,stderr=subprocess.STDOUT,check=True)
print('CONFIGURE_BATCH_CHECKS_COMPLETE',flush=True)
