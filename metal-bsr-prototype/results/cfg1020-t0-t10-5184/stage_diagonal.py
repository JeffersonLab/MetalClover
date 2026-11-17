"""Prepare cfg1020 P(t,t), t=1..10, without modifying cfg1000 workers."""
from pathlib import Path
import ast
import datetime
import hashlib
import json
import shutil
import xml.etree.ElementTree as ET

ROOT = Path('/Users/admin/Documents/ChatGPT/MG for Apple GPU')
RESULTS = ROOT / 'metal-bsr-prototype/results'
REPORT = Path(__file__).resolve().parent
OLD = RESULTS / 'perambulator-t6-t10-5184'
WORK = Path('/private/tmp/chroma-cfg1020-production')
RUNTIME = Path('/private/tmp/chroma-times6-10-production/runtime')


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(4 * 1024**2), b''):
            digest.update(block)
    return digest.hexdigest()


WORK.mkdir(exist_ok=False)
baseline = RESULTS / 'perambulator-t10-5184'
manifest = json.loads((baseline / 'frozen-files.json').read_text())
for name, expected in manifest.items():
    assert sha(name) == expected, name
common = {name: expected for name, expected in manifest.items()
          if name.startswith(str(RUNTIME) + '/')}
gauge = json.loads((REPORT / 'gauge-metadata.json').read_text())
assert gauge['configuration_from_header'] == 1020
assert gauge['full_SciDAC_checksum_verified'] and gauge['FermionAction_identical_to_cfg1000']
assert sha(gauge['path']) == gauge['SHA256']
inputs = {item['time']: item for item in json.loads((REPORT / 'eigenvector-inputs.json').read_text())['inputs']}

driver = (OLD / 'batch_queue.py').read_text()
driver = driver.replace('cfg1000', 'cfg1020').replace('sdb1000', 'sdb1020')
driver = driver.replace("f'perambulator-t{t}-5184'", "f'cfg1020-perambulator-t{t}-5184'")
driver = driver.replace('assert 0 <= t <= 10', 'assert 1 <= t <= 10')
driver = driver.replace('assert 6 <= owner <= 10 and all(0 <= t <= owner',
                        'assert 1 <= owner <= 10 and all(1 <= t <= owner')
driver = driver.replace('assert 6 <= t <= 10', 'assert 1 <= t <= 10')
driver = driver.replace('choices=range(6, 11)', 'choices=range(1, 11)')
driver = driver.replace('chroma-times6-10-production/batch-storage-ready.json',
                        'chroma-cfg1020-production/batch-storage-ready.json')
driver = driver.replace("f'/private/tmp/chroma-time{t}-production'",
                        "f'/private/tmp/chroma-cfg1020-time{t}-production'")
driver = driver.replace('previous, _, _ = paths(t - 1)',
                        "previous = BASE / 'cfg1000-t10-unsmeared' if t == 1 else paths(t - 1)[0]")
driver = driver.replace("'physical_rhs': RHS, 'queue_PID'", "'configuration': 1020, 'physical_rhs': RHS, 'queue_PID'")
driver = driver.replace("waiting_for=f'time-{t-1} complete validation and finalization',",
                        "waiting_for=('cfg1000 time10' if t == 1 else f'cfg1020 time{t-1}') + ' complete validation/finalization and disk budget',")
start = driver.index('        _, projected = ensure_evidence(t - 1)')
stop = driver.index('        while shutil.disk_usage(BASE).free', start)
driver = driver[:start] + '''        if t == 1:
            old_report = RESULTS / 'perambulator-t10-5184'
            accepted = read_json(str(previous) + '.complete.json')
            proof = read_json(old_report / 'production-validation.json')
            old_standard = BASE / 'szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_10.sdb1000'
            projected = read_json(str(old_standard) + '.complete.json')
            assert accepted['status'] == projected['status'] == 'complete'
            assert accepted['source_time'] == accepted['sink_time'] == 10
            assert accepted['completed_rhs'] == accepted['physical_rhs'] == RHS
            assert accepted['session_count'] == accepted['MG_hierarchy_builds_across_sessions'] == 1
            assert bounded(accepted['maximum_full_relative']) and bounded(accepted['maximum_qprop_relative'])
            assert proof['status'] == 'pass' and proof['all_standard_record_values_readback_bitwise']
            assert proof['independent_DataLib_full_file_key_check'] == 'pass'
            assert projected['spin_records'] == 16 and bounded(projected['full_gamma5_relative_Frobenius'])
            assert old_standard.stat().st_size == projected['output_bytes']
            assert sha(old_standard) == projected['output_SHA256'] == proof['standard_perambulator_SHA256']
            write_json(report / 'cfg1000-time10-dependency.json',
                       {'status': 'pass', 'standard_SHA256': projected['output_SHA256'],
                        'all_original_full_FP64_checks_and_independent_reader_passed': True})
        else:
            _, projected = ensure_evidence(t - 1)
            shutil.copyfile(Path(str(previous) + '.queue.json'),
                            paths(t - 1)[2] / 'queue-complete.json')
        retired = retire_intermediates(t, range(1, t), 'dependency-retirement.json')
        for source_time in sorted(set(retired + ([t - 1] if t > 1 else []))):
            execute(f'render-dependency-time{source_time}',
                    ['python3', work / 'render_completed.py', source_time])
''' + driver[stop:]
driver = driver.replace("f'../perambulator-t{owner}-5184/{record_name}'", "f'../cfg1020-perambulator-t{owner}-5184/{record_name}'")
driver = driver.replace("assert proof.get('status') == 'pass'", "assert proof.get('status') == 'pass' and proof.get('configuration') == 1020")
driver = driver.replace("'owner_time': owner,", "'owner_time': owner, 'configuration': 1020,")
ast.parse(driver)
for path in [REPORT / 'batch_queue.py', WORK / 'batch_queue.py']:
    path.write_text(driver)
common[str(WORK / 'batch_queue.py')] = sha(WORK / 'batch_queue.py')
shutil.copyfile(REPORT / 'storage_watch.py', WORK / 'storage_watch.py')
common[str(WORK / 'storage_watch.py')] = sha(WORK / 'storage_watch.py')

capture = (OLD / 'capture_completed.py').read_text()
capture = capture.replace('cfg1000', 'cfg1020').replace('sdb1000', 'sdb1020')
capture = capture.replace('results/perambulator-t{t}-5184', 'results/cfg1020-perambulator-t{t}-5184')
capture = capture.replace('chroma-time{t}-production', 'chroma-cfg1020-time{t}-production')
start = capture.index('frozen_path=')
stop = capture.index('\nfrozen=json.loads', start)
capture = capture[:start] + "frozen_path=Path(f'/private/tmp/chroma-cfg1020-time{t}-production/frozen-files.json')" + capture[stop:]
capture = capture.replace("record={'status':'pass',", "record={'status':'pass','configuration':1020,")
(REPORT / 'capture_completed.py').write_text(capture)
renderer = (OLD / 'render_completed.py').read_text()
renderer = renderer.replace('Configuration 1000', 'Configuration 1020')
renderer = renderer.replace('results/perambulator-t{t}-5184', 'results/cfg1020-perambulator-t{t}-5184')
renderer = renderer.replace('chroma-time{t}-production', 'chroma-cfg1020-time{t}-production')
renderer = renderer.replace("staged = json.loads((local_staged if t>=6 else report / 'staged-run.json').read_text())",
                            'staged = json.loads(local_staged.read_text())')
start = renderer.index('if t == 2:')
stop = renderer.index("if validation.get('point_intermediate_retired", start)
renderer = renderer[:start] + renderer[stop:]
(REPORT / 'render_completed.py').write_text(renderer)

template = Path('/private/tmp/chroma-time10-production/full_run.py').read_text()
action_template = (baseline / 'action.xml').read_text()
base_plan = json.loads((baseline / 'staged-run.json').read_text())
jobs = []
for t in range(1, 11):
    item = inputs[t]
    local = Path(f'/private/tmp/chroma-cfg1020-time{t}-production')
    report = RESULTS / f'cfg1020-perambulator-t{t}-5184'
    local.mkdir(exist_ok=False)
    report.mkdir(exist_ok=False)
    (report / 'production').mkdir()
    action = action_template.replace('<t_sources>10</t_sources>', f'<t_sources>{t}</t_sources>')
    assert ET.tostring(ET.fromstring(action).find('.//FermionAction')) == ET.tostring(ET.fromstring(action_template).find('.//FermionAction'))
    (local / 'action.xml').write_text(action)
    shutil.copyfile(local / 'action.xml', report / 'action.xml')
    pipeline = template
    for before, after in [
        ('chroma-time10-production', f'chroma-cfg1020-time{t}-production'),
        ('cfg1000-t10-unsmeared', f'cfg1020-t{t}-unsmeared'),
        ('light.t0_10.sdb1000', f'light.t0_{t}.sdb1020'),
        ('eigs_mod.n5184/1000/eigs.t_10.mod1000', f'eigs_mod.n5184/1020/eigs.t_{t}.mod1020'),
        ('/Users/admin/Downloads/download/prop_tests/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265_per_cfg_1000.lime', gauge['path']),
        ("'--time','10'", f"'--time','{t}'"),
        ("'KEYS_ONLY','10'", f"'KEYS_ONLY','{t}'"),
        ("'source_time':10", f"'source_time':{t}"),
        ("'sink_time':10", f"'sink_time':{t}"),
        ("completed['source_time']==completed['sink_time']==10", f"completed['source_time']==completed['sink_time']=={t}"),
    ]:
        assert before in pipeline, before
        pipeline = pipeline.replace(before, after)
    assert 'cfg1000' not in pipeline and 'mod1000' not in pipeline
    pipeline = pipeline.replace("'status':'running','source_time'", "'status':'running','configuration':1020,'source_time'")
    ast.parse(pipeline)
    (local / 'full_run.py').write_text(pipeline)
    shutil.copyfile(local / 'full_run.py', report / 'production/full_run.py')
    frozen = dict(common)
    for name in ['capture_completed.py', 'render_completed.py']:
        shutil.copyfile(REPORT / name, local / name)
        shutil.copyfile(REPORT / name, report / 'production' / name)
        frozen[str(local / name)] = sha(local / name)
    for name in ['action.xml', 'full_run.py']:
        frozen[str(local / name)] = sha(local / name)
    plan = dict(base_plan)
    plan.update(configuration=1020, source_time=t, sink_time=t,
                dependency='cfg1000 time10 complete' if t == 1 else f'cfg1020 time{t-1} complete',
                gauge_file=gauge['path'], gauge_SHA256=gauge['SHA256'],
                eigenvector_file=item['path'], eigenvector_bytes=item['bytes'],
                eigenvector_SHA256=item['SHA256'], eigenvalue_count=5184,
                eigenvalue_order='descending, preserved exactly',
                eigenvector_header=dict(base_plan['eigenvector_header'], t_slice=t),
                runtime_action_XML=str(local / 'action.xml'), action_XML_SHA256=sha(local / 'action.xml'),
                action_change='cfg1020 header action verified identical; contraction source time changed',
                point_prefix=f'/Users/admin/Downloads/MetalClover-results/cfg1020-t{t}-unsmeared',
                standard_perambulator=f'/Users/admin/Downloads/MetalClover-results/szscl3_12_96_b1p50_t_x4p300_um0p0743_n1p265.prop.n5184.light.t0_{t}.sdb1020',
                staged_UTC=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                batch_storage_gate=str(WORK / 'batch-storage-ready.json'))
    for directory in [local, report]:
        (directory / 'frozen-files.json').write_text(json.dumps(frozen, indent=2) + '\n')
        (directory / 'staged-run.json').write_text(json.dumps(plan, indent=2) + '\n')
    jobs.append(plan)
(REPORT / 'jobs.json').write_text(json.dumps(jobs, indent=2) + '\n')
print('Staged cfg1020 P(t,t), t=1..10, using the unchanged private numerical runtime.')
