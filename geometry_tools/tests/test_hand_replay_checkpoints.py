"""Local integration guard: saved replay prefixes are never installable."""
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace
import unittest


ROOT=Path(__file__).resolve().parents[2]
BASE=ROOT/'.local-models/geometry-workbench/biodigital-male-circulation-v1'


class HandReplayCheckpointTests(unittest.TestCase):
    def test_incomplete_prefix_is_rejected_before_packing(self):
        first=BASE/'hand-right-guarded-v1566-s1/report.json'
        if not first.exists():self.skipTest('Local reviewed hand evidence not installed')
        report=json.loads(first.read_text())
        self.assertTrue(report['stageComplete'])
        self.assertFalse(report['complete'])
        self.assertEqual(report['rows'],[])
        spec=importlib.util.spec_from_file_location('compact_vessel_installer',ROOT/'scripts/install-compact-vessel-preview.py')
        module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
        output='test-forbidden-incomplete-hand-pack'
        self.assertFalse((BASE/output).exists())
        with self.assertRaises(AssertionError):
            module.prepare(SimpleNamespace(candidate=str(first),previous='blood-installed-v1560',output=output))
        self.assertFalse((BASE/output).exists())

    def test_final_replay_retains_every_saved_prefix_check(self):
        path=BASE/'hand-right-guarded-v1566-s4/report.json'
        if not path.exists():self.skipTest('Local reviewed hand evidence not installed')
        final=json.loads(path.read_text());self.assertTrue(final['complete'])
        self.assertEqual(len(final['trials']),24)
        for number in [1,2,3]:
            prefix=json.loads((BASE/f'hand-right-guarded-v1566-s{number}/report.json').read_text())
            self.assertEqual(final['trials'][:len(prefix['trials'])],prefix['trials'])
            self.assertEqual(final['fields'][:len(prefix['fields'])],prefix['fields'])
            self.assertEqual(final['sourceAtlasSha256'],prefix['sourceAtlasSha256'])
