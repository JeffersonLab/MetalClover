"""Synthetic unit fixtures for queue gates and destructive-action guards only."""

from pathlib import Path
import hashlib
import importlib.util
import json
import subprocess
import tempfile
import unittest
from unittest.mock import patch


class QueueSafety(unittest.TestCase):
    def setUp(self):
        spec = importlib.util.spec_from_file_location(
            'batch_queue_test', Path(__file__).with_name('batch_queue.py'))
        self.queue = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(self.queue)
        self.temporary = tempfile.TemporaryDirectory(prefix='metalclover-queue-unit-')
        root = Path(self.temporary.name)
        self.queue.BASE = root / 'outputs'
        self.queue.RESULTS = root / 'reports'
        self.queue.BASE.mkdir()
        self.queue.RESULTS.mkdir()

    def tearDown(self):
        self.temporary.cleanup()

    def fixture(self):
        """Tiny buffers stand in for storage; these are not numerical results."""
        q = self.queue
        q.POINT_BYTES, q.STANDARD_BYTES = 4, 8
        q.NUMERICAL_PAYLOAD_BYTES = 8
        point, standard, report = q.paths(5)
        report.mkdir()
        q.paths(6)[2].mkdir()
        payload = Path(str(point) + '.npy')
        payload.write_bytes(b'data')
        standard.write_bytes(b'retained')
        point_hash = hashlib.sha256(b'data').hexdigest()
        standard_hash = hashlib.sha256(b'retained').hexdigest()
        accepted = {'status': 'complete', 'source_time': 5, 'sink_time': 5,
                    'physical_rhs': 20736, 'completed_rhs': 20736,
                    'num_vectors': 5184, 'source_spins': 4, 'session_count': 1,
                    'MG_hierarchy_builds_across_sessions': 1,
                    'maximum_full_relative': 2e-14, 'maximum_qprop_relative': 2e-14,
                    'full_CPU_FP64_tolerance': 1e-12, 'payload_bytes': 4,
                    'payload_sha256': point_hash}
        projected = {'status': 'complete', 'source_time': 5, 'sink_time': 5,
                     'spin_records': 16, 'num_vectors': 5184,
                     'full_gamma5_relative_Frobenius': 1e-14,
                     'max_QDP_dot_scaled_error': 1e-15, 'output_bytes': 8,
                     'output_SHA256': standard_hash}
        proof = {'status': 'pass', 'configuration': 1020, 'source_time': 5, 'sink_time': 5,
                 'physical_rhs': 20736, 'source_vectors': 5184,
                 'hierarchies': 1, 'sessions': 1, 'native_array_calls': 1296,
                 'RHS_per_call': 16, 'point_SHA256': point_hash,
                 'standard_perambulator_SHA256': standard_hash,
                 'all_point_column_checksums_verified': True,
                 'all_standard_record_values_readback_bitwise': True,
                 'independent_DataLib_full_file_key_check': 'pass',
                 'unchanged_frozen_source_executable_library_hashes': True}
        q.write_json(str(point) + '.complete.json', accepted)
        q.write_json(str(standard) + '.complete.json', projected)
        q.write_json(str(point) + '.pipeline.json',
                     {'status': 'complete', 'source_time': 5, 'sink_time': 5})
        q.write_json(report / 'production-validation.json', proof)
        return point, payload, standard, accepted

    def test_missing_future_dependency_waits(self):
        q = self.queue
        self.assertEqual(q.read_json(q.BASE / 'not-started.json'), {})
        self.assertFalse(q.dependency_ready({}, {}))
        self.assertFalse(q.dependency_ready({'status': 'complete'}, {'status': 'queued'}))
        self.assertFalse(q.dependency_ready({'status': 'failed'}, {'status': 'complete'}))
        self.assertTrue(q.dependency_ready({'status': 'complete'}, {'status': 'complete'}))

    def test_residual_gate_rejects_nonfinite_and_large_values(self):
        q = self.queue
        for value in [float('nan'), float('inf'), -1e-20, 1.01e-12]:
            self.assertFalse(q.bounded(value))
        self.assertTrue(q.bounded(2.5e-14))

    def test_unattended_storage_plan_must_cover_entire_batch(self):
        q = self.queue
        self.assertFalse(q.storage_ready({}))
        good = {'status': 'pass', 'forecast_all_runs_fit': True,
                'all_standard_perambulators_retained': True,
                'minimum_forecast_free_GiB': 7.0,
                'minimum_new_job_start_free_GiB': 20.0}
        self.assertTrue(q.storage_ready(good))
        self.assertFalse(q.storage_ready(dict(good, minimum_new_job_start_free_GiB=17.0)))
        self.assertFalse(q.storage_ready(dict(good, minimum_forecast_free_GiB=3.9)))

    def test_incomplete_solve_cannot_retire(self):
        point, payload, standard, accepted = self.fixture()
        accepted['status'] = 'running'
        self.queue.write_json(str(point) + '.complete.json', accepted)
        with self.assertRaises(AssertionError):
            self.queue.retire_intermediates(6, [5], 'retirement.json')
        self.assertEqual(payload.read_bytes(), b'data')
        self.assertEqual(standard.read_bytes(), b'retained')

    def test_wrong_time_cannot_retire(self):
        point, payload, standard, accepted = self.fixture()
        accepted['source_time'] = 4
        self.queue.write_json(str(point) + '.complete.json', accepted)
        with self.assertRaises(AssertionError):
            self.queue.retire_intermediates(6, [5], 'retirement.json')
        self.assertEqual(payload.read_bytes(), b'data')
        self.assertEqual(standard.read_bytes(), b'retained')

    def test_changed_payload_cannot_retire(self):
        _, payload, standard, _ = self.fixture()
        payload.write_bytes(b'edit')
        with self.assertRaises(AssertionError):
            self.queue.retire_intermediates(6, [5], 'retirement.json')
        self.assertEqual(payload.read_bytes(), b'edit')
        self.assertEqual(standard.read_bytes(), b'retained')

    def test_only_verified_intermediate_is_removed(self):
        _, payload, standard, _ = self.fixture()
        closed = subprocess.CompletedProcess([], 1, '', '')
        with patch.object(self.queue.subprocess, 'run', return_value=closed):
            self.assertEqual(self.queue.retire_intermediates(6, [5], 'retirement.json'), [5])
        self.assertFalse(payload.exists())
        self.assertEqual(standard.read_bytes(), b'retained')
        record = self.queue.read_json(self.queue.paths(6)[2] / 'retirement.json')
        self.assertEqual(record['status'], 'complete')
        self.assertTrue(record['files'][0]['retired'])

    def test_wrong_configuration_cannot_retire(self):
        _, payload, standard, _ = self.fixture()
        report = self.queue.paths(5)[2]
        proof = self.queue.read_json(report / 'production-validation.json')
        proof['configuration'] = 1000
        self.queue.write_json(report / 'production-validation.json', proof)
        with self.assertRaises(AssertionError):
            self.queue.retire_intermediates(6, [5], 'retirement.json')
        self.assertTrue(payload.exists() and standard.exists())

    def test_future_or_out_of_scope_retirement_is_rejected(self):
        with self.assertRaises(AssertionError):
            self.queue.retire_intermediates(6, [7], 'retirement.json')
        with self.assertRaises(AssertionError):
            self.queue.retire_intermediates(11, [5], 'retirement.json')


if __name__ == '__main__':
    unittest.main()
