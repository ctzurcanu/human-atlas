"""Actual whole-surface movement closure and unrelated-source isolation."""
import importlib.util,json,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('forest',ROOT/'scripts/build-spatial-dependency-forest.py')
forest=importlib.util.module_from_spec(spec);spec.loader.exec_module(forest)

class PancreaticClosure(unittest.TestCase):
    def test_every_region_selection_carries_the_gland_and_its_ducts(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_bytes())
        f=forest.build_forest('male',a);ids=[p['id'] for p in a['parts']]
        regions=set(next(p for p in a['parts'] if p['id']=='ZA:Pancreas')['replacementParts'])
        expected=regions|set(a['pancreaticDuctInstallation']['partIds'])
        for selected in sorted(expected)+['ZA:Pancreas']:
            plan=forest.movement_plan(f,[selected],ids,np.empty((0,2),dtype=np.uint32))
            self.assertEqual(set(plan['movingPartIds']),expected)
            self.assertEqual(plan['mustMoveScaleTogetherGroups'],['assembly:pancreatic-exterior'])
            self.assertIn('ZA:Pancreatic duct',plan['movingPartIds'])
            self.assertIn('ZA:Accessory pancreatic duct',plan['movingPartIds'])
            self.assertNotIn('ZA:Superior mesenteric vein',plan['movingPartIds'])
            self.assertNotIn('ZA:Pancreas',plan['movingPartIds'])
            self.assertNotIn('LOCAL:male:BOWELWALL:duodenum:serosa',plan['movingPartIds'])
    def test_shared_gland_gut_interfaces_are_explicit_refit_obligations(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json').read_bytes())
        f=forest.build_forest('male',a);ids=[p['id']for p in a['parts']]
        marker=a['finitePancreaticGlandIntegration'];self.assertEqual(len(marker['partIds']),20)
        self.assertEqual(sum(pid.startswith('LOCAL:male:BOWELWALL:duodenum:')for pid in marker['partIds']),12)
        plan=forest.movement_plan(f,['TA98PART:ZA:Pancreas:head'],ids,np.empty((0,2),dtype=np.uint32))
        interfaces=[row for row in f['interfaces']if row['kind']=='shared-tissue-interface']
        self.assertTrue(any({row['from'],row['to']}=={'TA98PART:ZA:Pancreas:head','LOCAL:male:BOWELWALL:duodenum:serosa'}for row in interfaces))
        self.assertTrue(any({row['from'],row['to']}=={'TA98PART:ZA:Pancreas:head','TA98PART:FJ3079:A05.8.02.013'}for row in interfaces))
        self.assertFalse(plan['geometryModified']);self.assertFalse(plan['cameraMovementAuthorized'])
        altered=json.loads(json.dumps(a));p=next(p for p in altered['parts']if p['id']=='LOCAL:male:BOWELWALL:duodenum:serosa');p['geometryOrigin']=[0,0,0]
        with self.assertRaises(AssertionError):forest.build_forest('male',altered)
    def test_a_cap_assembly_tag_alone_never_creates_mechanical_ownership(self):
        a={'parts':[dict(id='A',name='A',sourceId='A',sectionAssembly='TA98:pancreas'),
                    dict(id='B',name='B',sourceId='B',sectionAssembly='TA98:pancreas')]}
        f=forest.build_forest('male',a)
        self.assertNotEqual(f['partToGroup']['A'],f['partToGroup']['B'])
    def test_female_source_and_hidden_donors_remain_separate(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json').read_bytes())
        f=forest.build_forest('female',a);ids=[p['id'] for p in a['parts']]
        plan=forest.movement_plan(f,['LOCAL:female:pancreas:1039'],ids,np.empty((0,2),dtype=np.uint32))
        marker=a.get('finitePancreaticGlandIntegration')
        expected=set(marker['glandRegionIds'].values())|set(a['pancreaticDuctInstallation']['partIds'])if marker else{'LOCAL:female:pancreas:1039'}
        self.assertEqual(set(plan['movingPartIds']),expected)
        self.assertFalse(any(pid.startswith('HRA:') for pid in plan['movingPartIds']))
        if marker:
            self.assertEqual(len(marker['partIds']),38)
            self.assertNotIn('LOCAL:female:BOWELWALL:duodenum:serosa',plan['movingPartIds'])
            self.assertNotIn('LOCAL:female:cystic_duct:authored',plan['movingPartIds'])
            altered=json.loads(json.dumps(a));next(p for p in altered['parts']if p['id']=='LOCAL:female:bile_duct:1031')['geometryOrigin']=[0,0,0]
            with self.assertRaises(AssertionError):forest.build_forest('female',altered)

if __name__=='__main__':unittest.main()
