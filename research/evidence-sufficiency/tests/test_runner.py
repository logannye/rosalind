import copy
import json
from pathlib import Path
import sys
import unittest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from run_pilot import advancement_allowed, candidate_alts, locate, validate_protocol


class BoundaryTests(unittest.TestCase):
    def setUp(self):
        self.protocol = json.loads((HERE / 'protocol-p0.json').read_text())

    def test_candidates_are_truth_independent_and_include_each_eligible_alt(self):
        row = {'ref':'A','callable_depth':100,'a':94,'c':2,'g':3,'t':1}
        self.assertEqual(candidate_alts(row, self.protocol['candidate_rule']), ['C','G'])
        row['callable_depth'] = 101
        self.assertEqual(candidate_alts(row, self.protocol['candidate_rule']), ['G'])
        row['ref'] = 'N'
        self.assertEqual(candidate_alts(row, self.protocol['candidate_rule']), [])

    def test_roles_and_guard_coordinates(self):
        validate_protocol(self.protocol)
        intervals = self.protocol['intervals']
        self.assertEqual(locate(intervals,'chr20',30_494_999)['partition'],'calibration')
        with self.assertRaises(ValueError):
            locate(intervals,'chr20',30_495_000)
        self.assertEqual(locate(intervals,'chr20',30_505_000)['partition'],'threshold')
        bad = copy.deepcopy(self.protocol)
        bad['intervals'][3]['start'] = 30_494_999
        with self.assertRaisesRegex(ValueError,'overlapping'):
            validate_protocol(bad)

    def test_advancement_requires_recognized_development_gate_and_every_external_check(self):
        correct = dict.fromkeys(['independent_counts','native_parquet_features',
                                 'saved_only_replay','prediction_equivalence'],True)
        gate = {'status':'proceed_to_p1','proceed_to_p1':True,'test_outcomes_used':False}
        self.assertTrue(advancement_allowed(gate,correct))
        self.assertFalse(advancement_allowed({'passed':True},correct))
        for key in correct:
            self.assertFalse(advancement_allowed(gate,{**correct,key:False}))
        self.assertFalse(advancement_allowed({**gate,'test_outcomes_used':True},correct))
        self.assertFalse(advancement_allowed({**gate,'status':'insufficient_data'},correct))


if __name__ == '__main__':
    unittest.main()
