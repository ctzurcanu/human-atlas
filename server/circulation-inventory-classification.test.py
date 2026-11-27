"""Use actual authority rows to keep vessel recipients out of feature classes."""
import importlib.util
from pathlib import Path
import sqlite3
import unittest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('circulation_inventory', ROOT/'server/build-circulation-inventory.py')
inventory = importlib.util.module_from_spec(spec)
spec.loader.exec_module(inventory)


class VesselClassification(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with sqlite3.connect('file:'+str(inventory.AUTH)+'?mode=ro', uri=True) as conn:
            conn.row_factory = sqlite3.Row
            cls.terms = inventory.corrected_terms({r['id']: dict(r) for r in conn.execute('select * from ta98')})

    def test_named_recipient_and_landmark_do_not_turn_vessels_into_organs_or_curves(self):
        for code in ['A12.2.05.007', 'A12.2.05.008', 'A12.2.05.022', 'A12.2.08.014', 'A12.2.08.031', 'A12.3.04.005',
                     'A12.2.03.209', 'A12.2.06.009', 'A12.2.07.037', 'A12.2.07.064',
                     'A12.2.07.065', 'A12.2.07.066', 'A12.3.05.104', 'A12.3.11.020']:
            with self.subTest(code=code):
                self.assertEqual(inventory.classification(self.terms[code], {}), 'vessel-or-vascular-network')

    def test_wall_valve_and_root_spaces_keep_their_types(self):
        for code, expected in [('A12.0.00.017', 'wall-layer'), ('A12.0.00.018', 'wall-layer'),
                               ('A12.0.00.019', 'wall-layer'), ('A12.0.00.023', 'valve-or-leaflet'),
                               ('A12.2.03.002', 'space-or-anatomical-boundary')]:
            with self.subTest(code=code):
                self.assertEqual(inventory.classification(self.terms[code], {}), expected)

    def test_true_microvascular_structures_still_require_microvascular_geometry(self):
        rows = [r for r in self.terms.values() if r['name_en'].lower() in {'arteriole', 'venule', 'capillary', 'vasa vasorum'}]
        self.assertGreaterEqual(len(rows), 3)
        for row in rows:
            with self.subTest(code=row['id']):
                self.assertEqual(inventory.classification(row, {}), 'microscopic-tissue-or-network')


if __name__ == '__main__':
    unittest.main()
