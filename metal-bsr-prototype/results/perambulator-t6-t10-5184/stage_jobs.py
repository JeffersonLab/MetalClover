from pathlib import Path
import ast
import datetime
import hashlib
import json
import os
import shlex
import shutil
import subprocess
import xml.etree.ElementTree as ET

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
RESULTS = ROOT / 'metal-bsr-prototype/results'
BATCH = RESULTS / 'perambulator-t6-t10-5184'
WORK = Path('/private/tmp/chroma-times6-10-production')
RUNTIME = WORK / 'runtime'
BASELINE = RESULTS / 'perambulator-t5-5184'
baseline = json.loads((BASELINE / 'frozen-files.json').read_text())


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


for name, expected in baseline.items():
    if name.startswith('/'):
        assert sha(name) == expected, name
for directory in ['bin', 'libexec', 'lib', 'sources', 'controllers']:
    (RUNTIME / directory).mkdir(parents=True, exist_ok=True)
provenance = {}
frozen_common = {}


def snapshot(source, destination, role):
    source, destination = Path(source), Path(destination)
    shutil.copy2(source, destination)
    expected = sha(source)
    assert sha(destination) == expected
    provenance[str(source)] = {'snapshot': str(destination), 'SHA256': expected,
                               'bytes': source.stat().st_size, 'role': role}
    frozen_common[str(destination)] = expected


binary_dir = Path('/private/tmp/chroma-metal-validation/enabled/mainprogs/tests')
for name in ['t_metal_clover_unsmeared_time', 'project_perambulator_time', 'reader_check_time']:
    snapshot(binary_dir / name, RUNTIME / 'libexec' / name, 'unchanged native executable')
    wrapper = RUNTIME / 'bin' / name
    wrapper.write_text('#!/bin/sh\n'
                       'export DYLD_LIBRARY_PATH=' + shlex.quote(str(RUNTIME / 'lib')) + '\n'
                       'if [ "${METALCLOVER_DYLD_AUDIT:-0}" = "1" ]; then\n'
                       '  export DYLD_PRINT_LIBRARIES=1\nfi\n'
                       'exec ' + shlex.quote(str(RUNTIME / 'libexec' / name)) + ' "$@"\n')
    wrapper.chmod(0o755)
    frozen_common[str(wrapper)] = sha(wrapper)

libraries = [Path('/private/tmp/metalclover-install/lib/libmetalclover.dylib')]
visited = set()
to_scan = list((RUNTIME / 'libexec').iterdir())
while to_scan or libraries:
    if libraries:
        library = libraries.pop()
        if str(library) in visited:
            continue
        visited.add(str(library))
        target = RUNTIME / 'lib' / library.name
        if target.exists():
            assert sha(target) == sha(library), library
        else:
            snapshot(library, target, 'unchanged native shared library')
        to_scan.append(library)
    executable = to_scan.pop()
    dependencies = subprocess.check_output(['/usr/bin/otool', '-L', str(executable)], text=True)
    for line in dependencies.splitlines()[1:]:
        dependency = line.strip().split(' (', 1)[0]
        if dependency.startswith('@loader_path/'):
            dependency = str(executable.parent / dependency[len('@loader_path/'):])
        if dependency.startswith('/opt/homebrew/') and dependency not in visited:
            libraries.append(Path(dependency))

for source in [ROOT / 'metal-bsr-prototype/run_unsmeared.py',
               ROOT / 'metal-bsr-prototype/run_projected.py']:
    snapshot(source, RUNTIME / 'controllers' / source.name, 'unchanged controller')
for source in [Path('/Users/admin/Downloads/chroma.metalclover/mainprogs/tests/t_metal_clover_unsmeared.cc'),
               ROOT / 'metal-bsr-prototype/project_perambulator.cc',
               ROOT / 'metal-bsr-prototype/perambulator_reader_check.cc']:
    snapshot(source, RUNTIME / 'sources' / source.name, 'validated numerical source snapshot')

capture = Path('/private/tmp/chroma-time5-production/capture_completed.py').read_text()
old = "frozen=json.loads((dst/('validation/frozen-files.json' if t==1 else 'frozen-files.json')).read_text())"
new = "frozen_path=Path(f'/private/tmp/chroma-time{t}-production/frozen-files.json') if t>=6 else dst/('validation/frozen-files.json' if t==1 else 'frozen-files.json')\nfrozen=json.loads(frozen_path.read_text())"
assert old in capture
capture = capture.replace(old, new)
(WORK / 'capture_completed.py').write_text(capture)
shutil.copyfile(WORK / 'capture_completed.py', BATCH / 'capture_completed.py')

renderer = Path('/private/tmp/render_completed_peram.py').read_text()
old = "staged = json.loads((report / 'staged-run.json').read_text())"
new = "local_staged = Path(f'/private/tmp/chroma-time{t}-production/staged-run.json')\nstaged = json.loads((local_staged if t>=6 else report / 'staged-run.json').read_text())"
assert old in renderer
renderer = renderer.replace(old, new)
start = renderer.index("text += '''\nStorage: with explicit user approval,")
end = renderer.index("(report / 'README.md').write_text(text)", start)
renderer = renderer[:start] + '''if validation.get('point_intermediate_retired_after_full_validation'):
    retirement = validation['retirement_record']
    text += f"\\nThe accepted point intermediate was retired after complete verification and\\npoint/standard SHA256 rechecks. Its standard perambulator, all original inputs\\nand complete logs/residual evidence remain. See [the retirement record]({retirement}).\\n"
else:
    text += "\\nThe accepted point intermediate is retained pending the next validated handoff.\\n"
''' + renderer[end:]
(WORK / 'render_completed.py').write_text(renderer)
shutil.copyfile(WORK / 'render_completed.py', BATCH / 'render_completed.py')

shutil.copyfile(BATCH / 'batch_queue.py', WORK / 'batch_queue.py')
frozen_common[str(WORK / 'batch_queue.py')] = sha(WORK / 'batch_queue.py')

template = Path('/private/tmp/chroma-time5-production/full_run.py').read_text()
action_template = (BASELINE / 'action.xml').read_text()
action_fermion = ET.tostring(ET.fromstring(action_template).find('.//FermionAction'))
base_staged = json.loads((BASELINE / 'staged-run.json').read_text())
jobs = []
for t in range(6, 11):
    report = RESULTS / f'perambulator-t{t}-5184'
    work = Path(f'/private/tmp/chroma-time{t}-production')
    report.mkdir(exist_ok=False)
    work.mkdir(exist_ok=False)
    (report / 'production').mkdir()
    action = action_template.replace('<t_sources>5</t_sources>', f'<t_sources>{t}</t_sources>')
    assert ET.tostring(ET.fromstring(action).find('.//FermionAction')) == action_fermion
    (report / 'action.xml').write_text(action)
    shutil.copyfile(report / 'action.xml', work / 'action.xml')
    pipeline = template
    for before, after in [
        ('chroma-time5-production', f'chroma-time{t}-production'),
        ('cfg1000-t5-unsmeared', f'cfg1000-t{t}-unsmeared'),
        ('light.t0_5.sdb1000', f'light.t0_{t}.sdb1000'),
        ('eigs.t_5.mod1000', f'eigs.t_{t}.mod1000'),
        ('perambulator-t5-5184', f'perambulator-t{t}-5184'),
        ("'--time','5'", f"'--time','{t}'"),
        ("'KEYS_ONLY','5'", f"'KEYS_ONLY','{t}'"),
        ("'source_time':5", f"'source_time':{t}"),
        ("'sink_time':5", f"'sink_time':{t}"),
        ("completed['source_time']==completed['sink_time']==5",
         f"completed['source_time']==completed['sink_time']=={t}"),
    ]:
        assert before in pipeline, before
        pipeline = pipeline.replace(before, after)
    pipeline = pipeline.replace("bin=Path('/private/tmp/chroma-metal-validation/enabled/mainprogs/tests')",
                                f"bin=Path('{RUNTIME / 'bin'}')")
    pipeline = pipeline.replace("root/'metal-bsr-prototype/run_unsmeared.py'",
                                f"Path('{RUNTIME / 'controllers/run_unsmeared.py'}')")
    pipeline = pipeline.replace("root/'metal-bsr-prototype/run_projected.py'",
                                f"Path('{RUNTIME / 'controllers/run_projected.py'}')")
    pipeline = pipeline.replace(f"action=root/'metal-bsr-prototype/results/perambulator-t{t}-5184/action.xml'",
                                f"action=Path('{work / 'action.xml'}')")
    assert pipeline.count('==20736') == 1
    assert pipeline.count(f"completed['source_time']==completed['sink_time']=={t}") == 2
    assert 'chroma-metal-validation/enabled/mainprogs/tests' not in pipeline
    assert "root/'metal-bsr-prototype/run_" not in pipeline
    ast.parse(pipeline)
    (work / 'full_run.py').write_text(pipeline)
    shutil.copyfile(work / 'full_run.py', report / 'production/full_run.py')
    frozen = dict(frozen_common)
    for name in ['capture_completed.py', 'render_completed.py']:
        shutil.copyfile(WORK / name, work / name)
        shutil.copyfile(WORK / name, report / 'production' / name)
        frozen[str(work / name)] = sha(work / name)
    frozen[str(work / 'action.xml')] = sha(work / 'action.xml')
    frozen[str(work / 'full_run.py')] = sha(work / 'full_run.py')
    (report / 'frozen-files.json').write_text(json.dumps(frozen, indent=2) + '\n')
    shutil.copyfile(report / 'frozen-files.json', work / 'frozen-files.json')
    eig = Path(f'/Users/admin/Downloads/unitary_distillation/eigs_mod.n5184/1000/eigs.t_{t}.mod1000')
    with eig.open('rb') as stream:
        header = stream.read(1024)
    for token in [f'<t_slice>{t}</t_slice>'.encode(), b'<num_vecs>5184</num_vecs>',
                  b'<lattSize>12 12 12</lattSize>', b'<decay_dir>3</decay_dir>']:
        assert token in header, token
    staged = dict(base_staged)
    staged['eigenvector_header'] = dict(base_staged['eigenvector_header'], t_slice=t)
    staged.update(source_time=t, sink_time=t, dependency=f't={t-1} full pipeline and finalization',
                  eigenvector_file=str(eig), eigenvector_bytes=eig.stat().st_size,
                  eigenvector_SHA256=sha(eig), action_XML_SHA256=sha(work / 'action.xml'),
                  runtime_action_XML=str(work / 'action.xml'),
                  action_change=f'only Contractions/t_sources={t}; FermionAction unchanged',
                  point_prefix=f'/Users/admin/Downloads/MetalClover-results/cfg1000-t{t}-unsmeared',
                  standard_perambulator=f'/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb1000',
                  staged_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  isolated_runtime=str(RUNTIME),
                  minimum_free_disk_GiB_before_start=18.0,
                  batch_storage_gate=str(WORK / 'batch-storage-ready.json'),
                  intermediate_retention='capture all full checks, keep every standard perambulator, retire only accepted point intermediates with point/standard SHA256 rechecks',
                  free_disk_GiB_at_staging=shutil.disk_usage('/Users/admin/Downloads/MetalClover-results').free / 2**30)
    (report / 'staged-run.json').write_text(json.dumps(staged, indent=2) + '\n')
    shutil.copyfile(report / 'staged-run.json', work / 'staged-run.json')
    progress = Path('/private/tmp/chroma-time5-production/progress.py').read_text()
    progress = progress.replace("POINT = BASE / 'cfg1000-t5-unsmeared'", f"POINT = BASE / 'cfg1000-t{t}-unsmeared'")
    progress = progress.replace('light.t0_5.sdb1000', f'light.t0_{t}.sdb1000').replace("'time': 5", f"'time': {t}")
    (work / 'progress.py').write_text(progress)
    jobs.append(staged)
(BATCH / 'runtime-provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
(BATCH / 'jobs.json').write_text(json.dumps(jobs, indent=2) + '\n')
print(json.dumps({'times': [job['source_time'] for job in jobs],
                  'physical_rhs_total': 5 * 20736,
                  'isolated_runtime': str(RUNTIME), 'runtime_files': len(frozen_common)}, indent=2))
