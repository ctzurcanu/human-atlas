import sys,unittest,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
from open3d_source_identity import source_identity,native_representation_exists,sided_import_needed,retain_unpartitioned_branch_source
class SourceIdentity(unittest.TestCase):
    def test_actual_native_coverage_blocks_publisher_spelling_duplicate(self):
        atlas=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_text())
        native=[p['sourceId']for p in atlas['parts']if p['id'].startswith('ZA:')and p.get('sourceId')]
        for side in ['r','l']:
            self.assertTrue(native_representation_exists('Musculocutaneus nerve.'+side,native))
    def test_right_donor_can_fill_missing_left_without_replacing_right(self):
        right={source_identity('Musculocutaneous nerve.r')}
        self.assertTrue(sided_import_needed('Musculocutaneus nerve.r',right))
        both=right|{source_identity('Musculocutaneous nerve.l')}
        self.assertFalse(sided_import_needed('Musculocutaneus nerve.r',both))
    def test_unpartitioned_motor_branch_source_is_not_discarded_by_name(self):
        self.assertTrue(retain_unpartitioned_branch_source('Musculocutaneus nerve.r'))
        self.assertTrue(retain_unpartitioned_branch_source('Musculocutaneus nerve.l'))
        self.assertFalse(retain_unpartitioned_branch_source('Median nerve.r'))
    def test_missing_side_stays_importable_and_composite_stays_distinct(self):
        existing=['Musculocutaneous nerve.l']
        self.assertFalse(native_representation_exists('Musculocutaneus nerve.r',existing))
        self.assertTrue(native_representation_exists('Musculocutaneus nerve.l',existing))
        self.assertFalse(native_representation_exists('Musculocutaneus nerve - lateral antebrachial cutaneous nerve.l',existing))
        self.assertEqual(source_identity('Median nerve.r'),'median nerve.r')
if __name__=='__main__':unittest.main()
