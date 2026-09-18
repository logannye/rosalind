import importlib.util
import json
from pathlib import Path
import tempfile
import unittest

spec = importlib.util.spec_from_file_location('cgroup_probe', Path(__file__).with_name('cgroup_probe.py'))
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)


class OutcomeTests(unittest.TestCase):
    def result(self, **values):
        result = dict(analysis_exit_code=0, docker_oom_killed=False, events_before={'oom_kill': 0},
                      events_after={'oom_kill': 0}, timed_out=False, memory_max_bytes=128 << 20,
                      swap_max_bytes=0, success_artifact=True, success_receipt=True,
                      partial_artifact=False, partial_receipt=False, receipt_status='completed',
                      verified=True, assurance='cgroup-v2')
        result.update(values)
        return result

    def check(self, expected, result, **values):
        return probe.check_case(dict(expected=expected, hard_limit_bytes=128 << 20, **values), result)

    def test_completed_requires_verified_receipt_and_real_limit(self):
        self.assertEqual(self.check('completed', self.result(), require_os_limit=True), [])
        for change in [dict(verified=False), dict(memory_max_bytes=256 << 20),
                       dict(swap_max_bytes=None), dict(assurance='observed-only')]:
            self.assertTrue(self.check('completed', self.result(**change), require_os_limit=True))

    def test_signal_exit_alone_does_not_establish_oom(self):
        result = self.result(analysis_exit_code=137, success_artifact=False, success_receipt=False)
        self.assertTrue(self.check('oom', result))
        result['events_after'] = {'oom_kill': 1}
        self.assertEqual(self.check('oom', result), [])
        result['events_after'] = None
        result['docker_oom_killed'] = True
        self.assertEqual(self.check('oom', result), [])

    def test_failed_runs_cannot_keep_success_destinations(self):
        for code in [3, 4, 137]:
            result = self.result(analysis_exit_code=code, docker_oom_killed=True)
            self.assertIn('failed native invocation left a successful destination', self.check('oom', result))

    def test_capacity_requires_identified_integrity_checked_partial(self):
        result = self.result(analysis_exit_code=4, success_artifact=False, success_receipt=False,
                             partial_artifact=True, partial_receipt=True,
                             partial_status='resource-failed', partial_integrity=True)
        self.assertEqual(self.check('capacity', result), [])
        result['partial_status'] = 'completed'
        self.assertTrue(self.check('capacity', result))
        result['partial_status'] = 'resource-failed'
        result['partial_integrity'] = False
        self.assertTrue(self.check('capacity', result))

    def test_refusal_is_distinct_from_runtime_breach_and_oom(self):
        result = self.result(analysis_exit_code=3, success_artifact=False, success_receipt=False)
        self.assertEqual(self.check('refused', result), [])
        result['analysis_exit_code'] = 4
        self.assertTrue(self.check('refused', result))

    def test_cram_capacity_preflight_requires_specific_error_and_no_output(self):
        result = self.result(analysis_exit_code=4, success_artifact=False,
                             success_receipt=False, capacity_preflight_reported=True)
        self.assertEqual(self.check('capacity-preflight', result), [])
        for change in [dict(capacity_preflight_reported=False), dict(partial_artifact=True),
                       dict(partial_receipt=True), dict(analysis_exit_code=3),
                       dict(docker_oom_killed=True)]:
            self.assertTrue(self.check('capacity-preflight', dict(result, **change)))

    def test_cram_capacity_preflight_error_is_retained_and_classified(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, value in [('analysis.exit_code', '4'), ('memory.max.before', str(128 << 20)),
                                ('memory.swap.max.before', '0'), ('memory.events.before', 'oom_kill 0'),
                                ('memory.events.after', 'oom_kill 0'),
                                ('analysis.stderr', 'error: CRAM container mean read length exceeds declared max_read_len')]:
                (root / name).write_text(value)
            case = {'expected': 'capacity-preflight', 'hard_limit_bytes': 128 << 20}
            result = probe.collect(root, case, {'State': {'ExitCode': 4, 'OOMKilled': False}}, False)
            self.assertTrue(result['passed'], result['errors'])
            self.assertTrue(result['capacity_preflight_reported'])
            self.assertIn('analysis.stderr', result['files'])

    def test_missing_events_are_unknown_and_metadata_reads_are_capped(self):
        with tempfile.TemporaryDirectory() as root:
            path = Path(root) / 'metadata'
            self.assertIsNone(probe.events(path))
            path.write_text('{"status": "completed"}')
            with self.assertRaises(ValueError):
                probe.read_json(path, limit=4)

    def test_first_party_failed_receipt_is_classified_by_status(self):
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            for name, value in [('analysis.exit_code', '4'), ('memory.max.before', str(128 << 20)),
                                ('memory.swap.max.before', '0'), ('memory.events.before', 'oom_kill 0'),
                                ('memory.events.after', 'oom_kill 0')]:
                (root / name).write_text(value)
            (root / 'result.tsv.partial').write_text('#partial\n')
            receipt = root / 'result.tsv.manifest.json'
            receipt.write_text(json.dumps({'params': {'run_status': 'resource-failed', 'manifest_blake3': 'test'}}))
            Path(str(receipt) + '.verify.json').write_text(json.dumps({
                'ok': False, 'trust': {'receipt_integrity': {'state': 'satisfied'}}}))
            case = {'expected': 'capacity', 'hard_limit_bytes': 128 << 20}
            result = probe.collect(root, case, {'State': {'ExitCode': 4, 'OOMKilled': False}}, False)
            self.assertTrue(result['passed'], result['errors'])
            self.assertFalse(result['success_receipt'])
            self.assertTrue(result['partial_receipt'])


if __name__ == '__main__':
    unittest.main()
