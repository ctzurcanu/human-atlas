"""Consolidation preserves legacy objects, entry points and job provenance."""
from contextlib import redirect_stdout
from importlib import import_module, util
from io import StringIO
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

from geometry_tools.common import ROOT, digest_file
from geometry_tools.jobs import implementation_signature


class ConsolidationTests(unittest.TestCase):
    def test_all_legacy_modules_are_aliases_of_one_maintained_implementation(self):
        old_paths = [p for p in (ROOT / 'geometry').glob('*.py') if p.stem != '__init__']
        self.assertTrue(old_paths)
        for old in old_paths:
            canonical = import_module('geometry_tools.algorithms.' + old.stem)
            legacy = import_module('geometry.' + old.stem)
            self.assertIs(legacy, canonical, old.stem)
            self.assertEqual(Path(legacy.__file__).resolve(), ROOT / 'geometry_tools/algorithms' / old.name)

    def test_bare_old_imports_work_outside_the_project_directory(self):
        with tempfile.TemporaryDirectory() as td:
            env = dict(os.environ, PYTHONPATH=str(ROOT / 'geometry'))
            result = subprocess.run([sys.executable, '-B', '-c',
                'from tube_sections import TubeGeometryError as old; '
                'from geometry_tools.algorithms.tube_sections import TubeGeometryError as new; '
                'import branch_junction_envelopes; assert old is new'],
                cwd=td, env=env, capture_output=True, text=True, timeout=30)
            self.assertEqual(result.returncode, 0, result.stderr)

    def test_explicit_file_loaders_and_private_helpers_keep_working(self):
        spec = util.spec_from_file_location('consolidation_legacy_cage', ROOT / 'geometry/longitudinal_organ_cage.py')
        module = util.module_from_spec(spec); spec.loader.exec_module(module)
        from geometry_tools.algorithms.longitudinal_organ_cage import LongitudinalOrganCage
        self.assertIs(module.LongitudinalOrganCage, LongitudinalOrganCage)
        from geometry.tube_tool import _full_mesh_checks as old
        from geometry_tools.algorithms.tube_tool import _full_mesh_checks as new
        self.assertIs(old, new)

    def test_both_cli_entry_points_remain_available(self):
        for module in ['tube_tool', 'audit_alimentary_atlas']:
            for package in ['geometry', 'geometry_tools.algorithms']:
                result = subprocess.run([sys.executable, '-B', '-m', package + '.' + module, '--help'],
                                        cwd=ROOT, capture_output=True, text=True, timeout=30)
                self.assertEqual(result.returncode, 0, result.stderr)
                self.assertIn('usage:', result.stdout)

    def test_alimentary_cli_default_root_survives_directory_move(self):
        from geometry_tools.algorithms import audit_alimentary_atlas as module
        report = dict(components=[], mouths=[], oesophagus_to_anus_pass=True, finite_anal_canal_lining_ids=[])
        with tempfile.TemporaryDirectory() as td:
            output = str(Path(td) / 'result.json')
            with patch.object(sys, 'argv', ['audit', '--catalogue', 'fixture.json', '--report', output]), \
                 patch.object(module, 'audit_catalogue', return_value=report) as run, redirect_stdout(StringIO()):
                self.assertEqual(module.main(), 0)
            self.assertEqual(run.call_args.args[1], str(ROOT))
            self.assertEqual(json.loads(Path(output).read_text()), report)

    def test_jobs_bind_nested_algorithm_sources_and_skip_test_files(self):
        algorithm = ROOT / 'geometry_tools/algorithms/tube_sections.py'
        seen = []
        def changed(path):
            seen.append(Path(path))
            return 'changed-algorithm' if Path(path) == algorithm else digest_file(path)
        before = implementation_signature()
        with patch('geometry_tools.jobs.digest_file', side_effect=changed):
            after = implementation_signature()
        self.assertIn(algorithm, seen); self.assertNotEqual(before, after)
        self.assertFalse(any('tests' in p.parts for p in seen))


if __name__ == '__main__':
    unittest.main()
