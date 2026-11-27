from pathlib import Path
import importlib.util,json,sys,unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('suprarenal_forest',ROOT/'scripts/build-spatial-dependency-forest.py');F=importlib.util.module_from_spec(spec);spec.loader.exec_module(F)

class SuprarenalClosure(unittest.TestCase):
    def test_finite_glands_and_retired_sources_move_together(self):
        for sex in ['female','male']:
            with self.subTest(sex=sex):
                a=json.loads((ROOT/f'.local-models/ta98-review/bowel-wall-local-atlas/{sex}/atlas.json').read_bytes());marker=a['finiteSuprarenalIntegration'];forest=F.build_forest(sex,a)
                self.assertEqual(len(marker['partIds']),6)
                for g in marker['glands']:
                    expected='assembly:suprarenal:'+g['side'];self.assertEqual({forest['partToGroup'][pid]for pid in g['partIds']+[g['sourcePartId']]},{expected})
                kidney='LOCAL:female:l_kidney:189'if sex=='female'else 'ZA:Kidney.l';self.assertNotEqual(forest['partToGroup'][kidney],'assembly:suprarenal:left')
    def test_reject_split_registered_bank(self):
        a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/female/atlas.json').read_bytes());pid=a['finiteSuprarenalIntegration']['partIds'][0];next(p for p in a['parts']if p['id']==pid)['geometryOrigin']=[0,0,0]
        with self.assertRaises(AssertionError):F.build_forest('female',a)
if __name__=='__main__':unittest.main()
