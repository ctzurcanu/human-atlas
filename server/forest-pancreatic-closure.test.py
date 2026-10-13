"""Actual whole-surface movement closure and unrelated-source isolation."""
import importlib.util,json,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('forest',ROOT/'scripts/build-spatial-dependency-forest.py')
forest=importlib.util.module_from_spec(spec);spec.loader.exec_module(forest)

class PancreaticClosure(unittest.TestCase):
    def test_every_native_region_selection_carries_the_whole_active_exterior(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_bytes())
        f=forest.build_forest('male',a);ids=[p['id'] for p in a['parts']]
        expected=set(next(p for p in a['parts'] if p['id']=='ZA:Pancreas')['replacementParts'])
        for selected in sorted(expected)+['ZA:Pancreas']:
            plan=forest.movement_plan(f,[selected],ids,np.empty((0,2),dtype=np.uint32))
            self.assertEqual(set(plan['movingPartIds']),expected)
            self.assertEqual(plan['mustMoveScaleTogetherGroups'],['assembly:pancreatic-exterior'])
            self.assertNotIn('ZA:Pancreatic duct',plan['movingPartIds'])
            self.assertNotIn('ZA:Accessory pancreatic duct',plan['movingPartIds'])
            self.assertNotIn('ZA:Superior mesenteric vein',plan['movingPartIds'])
            self.assertNotIn('ZA:Pancreas',plan['movingPartIds'])
    def test_a_cap_assembly_tag_alone_never_creates_mechanical_ownership(self):
        a={'parts':[dict(id='A',name='A',sourceId='A',sectionAssembly='TA98:pancreas'),
                    dict(id='B',name='B',sourceId='B',sectionAssembly='TA98:pancreas')]}
        f=forest.build_forest('male',a)
        self.assertNotEqual(f['partToGroup']['A'],f['partToGroup']['B'])
    def test_female_source_and_hidden_donors_remain_separate(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json').read_bytes())
        f=forest.build_forest('female',a);ids=[p['id'] for p in a['parts']]
        plan=forest.movement_plan(f,['LOCAL:female:pancreas:1039'],ids,np.empty((0,2),dtype=np.uint32))
        self.assertEqual(plan['movingPartIds'],['LOCAL:female:pancreas:1039'])
        self.assertFalse(any(pid.startswith('HRA:') for pid in plan['movingPartIds']))

if __name__=='__main__':unittest.main()
