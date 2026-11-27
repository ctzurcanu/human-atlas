"""Real local orbital joins: full boundaries, sidedness and separation guards."""
import importlib.util,json,unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.source_port_dependencies import check_source_mouth

B=Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')

class OcularMouthTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.release=B/'blood-ocular-source-mouths-v2715';prep=json.loads((cls.release/'prepare-report.json').read_text());cls.old=json.loads(Path(prep['sourceCatalogue']).read_text());cls.new=json.loads((cls.release/'circulation-atlas.json').read_text());cls.h=cls.new['bloodOcularSourceMouthHistory'][-1];cls.index=json.loads(Path(cls.h['indexPath']).read_text());cls.reader=AtlasReader();cls.parts={p['id']:p for p in cls.old['parts']};cls.meshes={}
        for r in cls.index['interfaces']:
            for k in ['first','second']:
                if r[k] not in cls.meshes:cls.meshes[r[k]]=cls.reader.read_part(cls.old,cls.parts[r[k]],frame='model:local-ta98:assembled-metres')
    def test_all_eight_complete_opposed_mouths_pass_without_looser_tolerances(self):
        self.assertEqual(len(self.index['interfaces']),8)
        for r in self.index['interfaces']:
            proof=check_source_mouth(self.meshes[r['first']],self.meshes[r['second']],r['firstRing'],r['secondRing'])
            self.assertTrue(proof['passed']);self.assertEqual(proof,r['proof']);self.assertLessEqual(proof['maximumGap_m'],1e-6)
    def test_installation_keeps_all_meshes_and_existing_ledger_literal(self):
        for key in ['parts','chunks','concepts','ta98TermBindings','materials','defaultHidden']:self.assertEqual(self.old[key],self.new[key])
        prep=json.loads((self.release/'prepare-report.json').read_text());old=json.loads((Path(prep['sourceCatalogue']).parent/'source-port-dependencies.json').read_text());new=json.loads((self.release/'source-port-dependencies.json').read_text());self.assertEqual(new['interfaces'][:len(old['interfaces'])],old['interfaces']);self.assertEqual(len(new['interfaces'])-len(old['interfaces']),8)
    def test_partial_contour_is_rejected(self):
        r=self.index['interfaces'][0]
        with self.assertRaises(ValueError):check_source_mouth(self.meshes[r['first']],self.meshes[r['second']],r['firstRing'][:-1],r['secondRing'])
    def test_opposite_eye_cannot_supply_a_same_side_join(self):
        r=self.index['interfaces'][0];other=r['second'].replace('-left_','-right_')
        proof=check_source_mouth(self.meshes[r['first']],self.meshes[other],r['firstRing'],r['secondRing']);self.assertFalse(proof['passed'])
    def test_separated_vessel_cannot_reuse_a_previously_passing_connection(self):
        r=self.index['interfaces'][0];m=self.meshes[r['first']];q=Mesh(m.vertices+np.array([0.,.0001,0.]),m.faces,m.frame,m.attributes,m.source)
        self.assertFalse(check_source_mouth(q,self.meshes[r['second']],r['firstRing'],r['secondRing'])['passed'])

if __name__=='__main__':unittest.main()
