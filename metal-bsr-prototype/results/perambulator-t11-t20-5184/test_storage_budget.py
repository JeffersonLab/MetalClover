"""Synthetic allocation forecasts; no numerical files or GPU work."""
import importlib.util
from pathlib import Path
import unittest

spec = importlib.util.spec_from_file_location('storage_watch', Path(__file__).with_name('storage_watch.py'))
storage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(storage)


class StorageBudget(unittest.TestCase):
    def states(self):
        return [{'time': t, 'queue_status': 'running' if t == 7 else 'queued',
                 'pipeline_exists': t == 7,
                 'point': [123, storage.POINT, storage.POINT] if t == 7 else [None, 0, 0],
                 'standard': [None, 0, 0]} for t in range(7, 21)]

    def test_current_capacity_does_not_enable_batch(self):
        plan = storage.forecast(35 * storage.GIB, self.states())
        self.assertFalse(plan['forecast_all_runs_fit'])
        self.assertGreater(plan['additional_free_space_required_GiB'], 60)

    def test_sufficient_capacity_preserves_every_standard(self):
        plan = storage.forecast(105 * storage.GIB, self.states())
        self.assertTrue(plan['forecast_all_runs_fit'])
        self.assertTrue(plan['all_standard_perambulators_retained'])
        self.assertGreaterEqual(plan['minimum_new_job_start_free_GiB'], 18)
        self.assertGreaterEqual(plan['minimum_forecast_free_GiB'], 4)

    def test_existing_allocations_are_not_charged_twice(self):
        states = self.states()
        first = storage.forecast(105 * storage.GIB, states)
        states[0]['standard'] = [456, storage.STANDARD, storage.STANDARD]
        second = storage.forecast(105 * storage.GIB - storage.STANDARD, states)
        self.assertAlmostEqual(first['final_forecast_free_GiB'], second['final_forecast_free_GiB'])

    def test_completed_outputs_stay_charged_to_measured_free_space(self):
        states = self.states()
        for state in states:
            state['queue_status'] = 'complete'
        plan = storage.forecast(20 * storage.GIB, states)
        self.assertTrue(plan['forecast_all_runs_fit'])
        self.assertEqual(plan['final_forecast_free_GiB'], 19)


if __name__ == '__main__':
    unittest.main()
