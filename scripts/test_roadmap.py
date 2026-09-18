"""Regression checks for roadmap validation before GitHub mutations."""
import copy
import json
import unittest
from unittest.mock import patch

import roadmap


class RoadmapTest(unittest.TestCase):
    def setUp(self):
        self.data = json.loads(roadmap.SOURCE.read_text())

    def test_repository_inventory_is_valid(self):
        roadmap.validate(self.data)

    def test_duplicate_mapping_is_rejected(self):
        self.data['tasks'][1]['issue_number'] = self.data['tasks'][0]['issue_number']
        self.data['tasks'][1]['issue_url'] = self.data['tasks'][0]['issue_url']
        with self.assertRaisesRegex(ValueError, 'duplicate roadmap issue'):
            roadmap.validate(self.data)

    def test_foreign_issue_url_is_rejected(self):
        self.data['tasks'][0]['issue_url'] = 'https://github.com/another/repo/issues/106'
        with self.assertRaisesRegex(ValueError, 'issue URL'):
            roadmap.validate(self.data)

    def test_completion_needs_evidence(self):
        self.data['tasks'][0].update(status='complete', evidence=[])
        with self.assertRaisesRegex(ValueError, 'completion needs evidence'):
            roadmap.validate(self.data)

    def test_conflicting_remote_mapping_refuses_before_writes(self):
        task = self.data['tasks'][0]
        remote = [{'number': task['issue_number'] + 1000, 'title': '[G01] Changed number'}]
        with patch.object(roadmap, 'gh_json', return_value=remote) as api, \
                patch.object(roadmap.subprocess, 'run') as subprocess:
            with self.assertRaisesRegex(ValueError, 'conflicting issue mapping'):
                roadmap.sync(self.data)
            self.assertEqual(api.call_count, 1)
            subprocess.assert_not_called()

    def test_sync_preserves_remote_labels_and_assignees(self):
        task = copy.deepcopy(self.data['tasks'][0])
        data = dict(self.data, tasks=[task], milestones={'G': self.data['milestones']['G']})
        calls = []
        def api(args, payload=None):
            calls.append((args, payload))
            if args[0] == 'issue':
                return [{'number': task['issue_number'], 'title': '[G01] Inventory'}]
            if 'milestones?' in args[1]:
                return [{'title': data['milestones']['G'], 'number': 1}]
            if payload is None:
                return {'labels': [{'name': 'maintainer-note'}], 'assignees': [{'login': 'collaborator'}]}
            return {'number': task['issue_number'], 'html_url': task['issue_url']}
        with patch.object(roadmap, 'gh_json', side_effect=api), \
                patch.object(roadmap.subprocess, 'run'), patch.object(roadmap, 'SOURCE'):
            roadmap.sync(data)
        payload = next(body for _, body in calls if body and 'labels' in body)
        self.assertIn('maintainer-note', payload['labels'])
        self.assertEqual(payload['assignees'], ['collaborator'])


if __name__ == '__main__':
    unittest.main()
