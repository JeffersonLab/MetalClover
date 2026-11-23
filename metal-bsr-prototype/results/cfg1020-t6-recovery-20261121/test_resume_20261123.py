"""Recovery and resource-control fixtures; never open production payloads."""
from pathlib import Path
import hashlib
import importlib.util
import tempfile
import unittest

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('resume', HERE / 'resume_20261123.py')
resume = importlib.util.module_from_spec(spec)
spec.loader.exec_module(resume)


class RecoverySafety(unittest.TestCase):
    def test_restore_is_exact_and_never_overwrites_present_mismatch(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'backup'
            target = Path(directory) / 'runtime' / 'executable'
            source.write_bytes(b'frozen numerical binary fixture')
            digest = hashlib.sha256(source.read_bytes()).hexdigest()
            self.assertTrue(resume.restore(source, target, digest, 0o755))
            self.assertEqual(target.read_bytes(), source.read_bytes())
            self.assertEqual(target.stat().st_mode & 0o777, 0o755)
            self.assertFalse(resume.restore(source, target, digest, 0o755))
            target.write_bytes(b'unexpected change')
            with self.assertRaises(AssertionError):
                resume.restore(source, target, digest, 0o755)
            self.assertEqual(target.read_bytes(), b'unexpected change')

    def test_corrupt_backup_cannot_be_published(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / 'backup'
            target = Path(directory) / 'runtime'
            source.write_bytes(b'corrupt')
            with self.assertRaises(AssertionError):
                resume.restore(source, target, '0' * 64)
            self.assertFalse(target.exists())

    def test_safe_prefix_runs_while_later_jobs_wait_for_space(self):
        storage = resume.module('storage_fixture', resume.BATCH / 'storage_watch.py')
        queue = resume.module('queue_fixture', resume.BATCH / 'batch_queue.py')
        states = [{'configuration': 1020, 'time': t,
                   'queue_status': 'failed' if t == 6 else 'queued',
                   'pipeline_exists': t == 6,
                   'point': [123, storage.POINT, storage.POINT] if t == 6 else [None, 0, 0],
                   'standard': [None, 0, 0]} for t in range(6, 11)]
        full = storage.forecast(30 * storage.GIB, states)
        self.assertFalse(full['forecast_all_runs_fit'])
        for t in range(7, 11):
            scoped = storage.forecast(30 * storage.GIB, [s for s in states if s['time'] <= t])
            self.assertEqual(queue.storage_ready(scoped), t <= 8)
            self.assertTrue(scoped['all_standard_perambulators_retained'])


if __name__ == '__main__':
    unittest.main()
