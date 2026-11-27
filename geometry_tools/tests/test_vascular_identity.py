"""Integration guards against ambiguous occurrences and stale source review."""
import copy
import json
import tempfile
import unittest
from pathlib import Path

from geometry_tools.vascular_identity import apply_reviewed_identities


class VascularIdentityTests(unittest.TestCase):
    def setUp(self):
        self.base = Path('.local-models/geometry-workbench/biodigital-male-circulation-v1').resolve()
        self.catalogue = self.base / 'blood-installed-v1455/circulation-atlas.json'
        self.source = json.loads(self.catalogue.read_text())
        self.plan_path = self.base / 'pancreatic-venous-identity-review-v1457/plan.json'
        self.plan = json.loads(self.plan_path.read_text())

    def test_same_name_with_another_recipient_cannot_be_mislabeled(self):
        plan = copy.deepcopy(self.plan)
        plan['bindings'][0]['term'] = 'A12.3.12.022'
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'plan.json'
            path.write_text(json.dumps(plan))
            with self.assertRaises(AssertionError):
                apply_reviewed_identities(copy.deepcopy(self.source), self.source, self.catalogue, path)

    def test_stale_render_review_cannot_be_installed(self):
        plan = copy.deepcopy(self.plan)
        plan['reviewSha256'] = '0' * 64
        with tempfile.TemporaryDirectory() as folder:
            path = Path(folder) / 'plan.json'
            path.write_text(json.dumps(plan))
            with self.assertRaises(AssertionError):
                apply_reviewed_identities(copy.deepcopy(self.source), self.source, self.catalogue, path)

    def test_identity_replay_changes_no_geometry_or_visibility(self):
        after = copy.deepcopy(self.source)
        result = apply_reviewed_identities(after, self.source, self.catalogue, self.plan_path)
        self.assertEqual(result['terms'], ['A12.3.12.029'])
        selected = set(result['partIds'])
        tags = {'ta98Term', 'ta98Kind', 'coverageLimitation'}
        for old, new in zip(self.source['parts'], after['parts']):
            if old['id'] in selected:
                self.assertEqual({k:v for k,v in old.items() if k not in tags},
                                 {k:v for k,v in new.items() if k not in tags})
            else:
                self.assertEqual(old, new)
        for key in ['chunks', 'materials', 'defaultHidden']:
            self.assertEqual(after[key], self.source[key])


if __name__ == '__main__':
    unittest.main()
