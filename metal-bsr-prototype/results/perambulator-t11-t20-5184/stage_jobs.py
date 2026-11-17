"""Stage a separate extension without modifying any existing worker/runtime."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import re
import shutil
import xml.etree.ElementTree as ET

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
RESULTS = ROOT / 'metal-bsr-prototype/results'
REPORT = Path(__file__).resolve().parent
WORK = Path('/private/tmp/chroma-times11-20-production')
OLD = RESULTS / 'perambulator-t6-t10-5184'
RUNTIME = Path('/private/tmp/chroma-times6-10-production/runtime')


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            h.update(block)
    return h.hexdigest()


WORK.mkdir(exist_ok=False)
inputs = json.loads((REPORT / 'eigenvector-inputs.json').read_text())['inputs']
baseline = RESULTS / 'perambulator-t10-5184'
frozen = json.loads((baseline / 'frozen-files.json').read_text())
for name, expected in frozen.items():
    assert sha(name) == expected, name
common = {name: value for name, value in frozen.items()
          if name.startswith(str(RUNTIME) + '/')}
assert len(common) >= 20
shutil.copyfile(REPORT / 'storage_watch.py', WORK / 'storage_watch.py')
common[str(WORK / 'storage_watch.py')] = sha(WORK / 'storage_watch.py')
driver = (OLD / 'batch_queue.py').read_text()
driver = driver.replace('assert 0 <= t <= 10', 'assert 0 <= t <= 20')
driver = driver.replace('assert 6 <= owner <= 10', 'assert 11 <= owner <= 20')
driver = driver.replace('assert 6 <= t <= 10', 'assert 11 <= t <= 20')
driver = driver.replace('choices=range(6, 11)', 'choices=range(11, 21)')
driver = driver.replace('chroma-times6-10-production/batch-storage-ready.json',
                        'chroma-times11-20-production/batch-storage-ready.json')
driver = driver.replace("waiting_for=f'time-{t-1} complete validation and finalization',",
                        "waiting_for=f'time-{t-1} complete validation/finalization and full-sequence disk budget',")
assert 'assert 11 <= t <= 20' in driver
ast.parse(driver)
for path in [REPORT / 'batch_queue.py', WORK / 'batch_queue.py']:
    path.write_text(driver)
common[str(WORK / 'batch_queue.py')] = sha(WORK / 'batch_queue.py')
for name in ['capture_completed.py', 'render_completed.py']:
    shutil.copyfile(OLD / name, REPORT / name)

template = Path('/private/tmp/chroma-time10-production/full_run.py').read_text()
action_template = (baseline / 'action.xml').read_text()
fermion_action = ET.tostring(ET.fromstring(action_template).find('.//FermionAction'))
base_plan = json.loads((baseline / 'staged-run.json').read_text())
jobs = []
for item in inputs:
    t = item['time']
    assert 11 <= t <= 20 and sha(item['path']) == item['SHA256']
    report = RESULTS / f'perambulator-t{t}-5184'
    local = Path(f'/private/tmp/chroma-time{t}-production')
    report.mkdir(exist_ok=False)
    local.mkdir(exist_ok=False)
    (report / 'production').mkdir()
    action = action_template.replace('<t_sources>10</t_sources>', f'<t_sources>{t}</t_sources>')
    assert ET.tostring(ET.fromstring(action).find('.//FermionAction')) == fermion_action
    (local / 'action.xml').write_text(action)
    shutil.copyfile(local / 'action.xml', report / 'action.xml')
    pipeline = template
    for before, after in [
        ('chroma-time10-production', f'chroma-time{t}-production'),
        ('cfg1000-t10-unsmeared', f'cfg1000-t{t}-unsmeared'),
        ('light.t0_10.sdb1000', f'light.t0_{t}.sdb1000'),
        ('eigs.t_10.mod1000', f'eigs.t_{t}.mod1000'),
        ("'--time','10'", f"'--time','{t}'"),
        ("'KEYS_ONLY','10'", f"'KEYS_ONLY','{t}'"),
        ("'source_time':10", f"'source_time':{t}"),
        ("'sink_time':10", f"'sink_time':{t}"),
        ("completed['source_time']==completed['sink_time']==10",
         f"completed['source_time']==completed['sink_time']=={t}"),
    ]:
        assert before in pipeline, before
        pipeline = pipeline.replace(before, after)
    assert pipeline.count(f"completed['source_time']==completed['sink_time']=={t}") == 2
    ast.parse(pipeline)
    (local / 'full_run.py').write_text(pipeline)
    shutil.copyfile(local / 'full_run.py', report / 'production/full_run.py')
    manifest = dict(common)
    for name in ['capture_completed.py', 'render_completed.py']:
        shutil.copyfile(REPORT / name, local / name)
        shutil.copyfile(REPORT / name, report / 'production' / name)
        manifest[str(local / name)] = sha(local / name)
    for name in ['action.xml', 'full_run.py']:
        manifest[str(local / name)] = sha(local / name)
    plan = dict(base_plan)
    plan.update(source_time=t, sink_time=t, dependency=f't={t-1} full pipeline and finalization',
                eigenvector_file=item['path'], eigenvector_bytes=item['bytes'],
                eigenvector_SHA256=item['SHA256'],
                eigenvector_header=dict(base_plan['eigenvector_header'], t_slice=t),
                runtime_action_XML=str(local / 'action.xml'),
                action_XML_SHA256=sha(local / 'action.xml'),
                action_change=f'only Contractions/t_sources={t}; FermionAction unchanged',
                point_prefix=f'/Users/admin/Downloads/MetalClover-results/cfg1000-t{t}-unsmeared',
                standard_perambulator=f'/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb1000',
                staged_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                batch_storage_gate=str(WORK / 'batch-storage-ready.json'))
    for directory in [report, local]:
        (directory / 'frozen-files.json').write_text(json.dumps(manifest, indent=2) + '\n')
        (directory / 'staged-run.json').write_text(json.dumps(plan, indent=2) + '\n')
    progress = Path('/private/tmp/chroma-time10-production/progress.py').read_text()
    progress = progress.replace('cfg1000-t10-unsmeared', f'cfg1000-t{t}-unsmeared')
    progress = progress.replace('light.t0_10.sdb1000', f'light.t0_{t}.sdb1000')
    progress = progress.replace("'time': 10", f"'time': {t}")
    (local / 'progress.py').write_text(progress)
    (report / 'README.md').write_text(f'''# Configuration 1000: source/sink t={t}

Queued extension: all 5,184 eigenvectors and four spins (20,736 RHS), after
complete time-{t-1} validation and finalization. Source and sink both equal {t}.
The action, numerical binaries and private runtime are unchanged. Caller/device
caps remain sixteen, one MG hierarchy is reused across the entire collection,
and both original full Chroma CPU FP64 checks remain at 1e-12 per RHS.

Input: `{item['path']}`, {item['bytes']:,} bytes, SHA256
`{item['SHA256']}`. The header reports 12³, decay_dir=3, t_slice={t},
num_vecs=5184. Supplied FP32 ordering/phases are preserved and promoted to FP64.

Standard output: `{plan['standard_perambulator']}`.
Sixteen 5,184×5,184 complex-FP64 records undergo complete readback, direct QDP
contractions, full gamma5 Hermiticity and the independent DataLib key reader.
Complete evidence is saved before retiring an accepted point intermediate;
point/standard SHA256 are rechecked. Every standard perambulator is retained.

The detached worker continues without a phone/chat connection, waiting for both
its predecessor and the storage guard. See [the batch record](../perambulator-t11-t20-5184/README.md)
for current disk-space status, manifests, scripts and the frozen runtime.
''')
    jobs.append(plan)
(REPORT / 'jobs.json').write_text(json.dumps(jobs, indent=2) + '\n')
print(json.dumps({'status': 'staged', 'times': [job['source_time'] for job in jobs],
                  'additional_RHS': 207360, 'unchanged_private_runtime': str(RUNTIME)}, indent=2))
