"""Current finite spleen owns its tissues, never the neighboring skeleton/vessels."""
import copy,importlib.util,json,sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'scripts'));spec=importlib.util.spec_from_file_location('spleen_forest',ROOT/'scripts/build-spatial-dependency-forest.py');F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)
class SpleenClosure(unittest.TestCase):
    def test_each_tissue_and_retired_surface_carry_one_complete_assembly(self):
        for sex in ['male','female']:
            a=json.loads((ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json').read_bytes());marker=a['finiteSpleenIntegration'];forest=F.build_forest(sex,a);expected=set(marker['partIds']);ids=[p['id']for p in a['parts']]
            for selected in [*expected,marker['sourcePartId']]:
                plan=F.movement_plan(forest,[selected],ids,np.empty((0,2),dtype=np.uint32));self.assertEqual(set(plan['movingPartIds']),expected);self.assertEqual(plan['mustMoveScaleTogetherGroups'],['assembly:spleen']);self.assertFalse(plan['cameraMovementAuthorized'])
    def test_incompatible_tissue_bank_is_rejected(self):
        for sex in ['male','female']:
            a=json.loads((ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json').read_bytes());pid=a['finiteSpleenIntegration']['partIds'][0];next(p for p in a['parts']if p['id']==pid)['geometryOrigin']=[0,0,0]
            with self.assertRaises(AssertionError):F.build_forest(sex,a)
if __name__=='__main__':unittest.main()
