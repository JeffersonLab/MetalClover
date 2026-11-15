"""Check private-library loading and actual argument lists without GPU work."""
from pathlib import Path
import ast
import hashlib
import json
import os
import py_compile
import subprocess

REPORT = Path(__file__).resolve().parent
WORK = Path('/private/tmp/chroma-times6-10-production')
RUNTIME = WORK / 'runtime'
checks = []
for name in ['t_metal_clover_unsmeared_time', 'project_perambulator_time']:
    wrapper = RUNTIME / 'bin' / name
    subprocess.run(['/bin/sh', '-n', str(wrapper)], check=True)
    environment = os.environ.copy()
    environment['METALCLOVER_DYLD_AUDIT'] = '1'
    result = subprocess.run([str(wrapper)], env=environment, capture_output=True, text=True)
    output = result.stdout + result.stderr
    (REPORT / (name + '.loader-audit.log')).write_text(output)
    assert result.returncode == 2 and 'Usage:' in output, output[-2000:]
    assert str(RUNTIME / 'lib/libmetalclover.dylib') in output
    assert 'NEAR_NULL_SVD' not in output
    checks.append({'name': name, 'exit_code': result.returncode,
                   'private_frozen_Metal_library_loaded': True,
                   'usage_only_no_operator_or_GPU_work': True})

for t in range(6, 11):
    work = Path(f'/private/tmp/chroma-time{t}-production')
    for name in ['full_run.py', 'capture_completed.py', 'render_completed.py', 'progress.py']:
        py_compile.compile(str(work / name), cfile=str(work / (name + 'c')), doraise=True)
    module = ast.parse((work / 'full_run.py').read_text())
    namespace = {}
    wanted = {'root', 'scratch', 'bin', 'point', 'peram', 'gauge', 'eigen', 'action', 'commands'}
    for node in module.body:
        if isinstance(node, (ast.Import, ast.ImportFrom)) or (
                isinstance(node, ast.Assign) and all(isinstance(target, ast.Name)
                                                   and target.id in wanted
                                                   for target in node.targets)):
            exec(compile(ast.Module(body=[node], type_ignores=[]), 'pipeline-commands', 'exec'), namespace)
    commands = dict(namespace['commands'])
    solve = list(map(str, commands['solve']))
    reader = list(map(str, commands['reader']))
    assert solve[solve.index('--time') + 1] == str(t)
    assert solve[solve.index('--vectors') + 1] == '5184'
    assert solve[solve.index('--max-rhs') + 1] == '16'
    assert reader[-1] == str(t)
    assert str(work / 'action.xml') == solve[solve.index('--action') + 1]
    assert f'eigs.t_{t}.mod1000' in solve[solve.index('--eigenvectors') + 1]
    frozen = json.loads((work / 'frozen-files.json').read_text())
    for name, expected in frozen.items():
        digest = hashlib.sha256()
        with Path(name).open('rb') as stream:
            for block in iter(lambda: stream.read(4 * 1024**2), b''):
                digest.update(block)
        assert digest.hexdigest() == expected, name
checks.append({'source_sink_times': list(range(6, 11)),
               'all_native_argument_lists_verified': True,
               'all_frozen_file_hashes_verified': True,
               'physical_RHS_per_run': 20736, 'caller_device_caps': 16})
(REPORT / 'runtime-validation.json').write_text(json.dumps({'status': 'pass', 'checks': checks}, indent=2) + '\n')
print(json.dumps(checks, indent=2))
