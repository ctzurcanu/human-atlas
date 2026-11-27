"""The anatomical heart assembly carries its real tissues, never the skeleton."""
from pathlib import Path
import importlib.util,json,sys,unittest
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
spec=importlib.util.spec_from_file_location('heart_forest',ROOT/'scripts/build-spatial-dependency-forest.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


class HeartForestTests(unittest.TestCase):
    def test_whole_heart_dependencies_in_both_current_models(self):
        for sex in ['male','female']:
            a=json.loads((ROOT/'.local-models/ta98-review/bowel-wall-local-atlas'/sex/'atlas.json').read_text())
            forest=module.build_forest(sex,a);m=a['finiteCardiacIntegration'];p={p['id']:p for p in a['parts']}
            g=next(g for root in forest['forest']for g in root['children']if g['id']=='assembly:heart')
            expected={pid for field in ['wallPartIds','valvePartIds','papillaryPartIds','chordalPartIds']for pid in m[field]}
            expected.update(pid for pid in m['rootPartIds']if p[pid].get('cardiacTemplateSourcePartId')or p[pid].get('cardiacRootWall'))
            self.assertEqual(set(g['activePartIds']),expected)
            self.assertTrue(g['moveScaleTogether'])
            self.assertFalse(any(p[pid]['system']=='skeletal'for pid in g['partIds']))
            if sex=='female':self.assertNotIn('LOCAL:female:aortic_arch:1513',g['partIds'])
            leaves=[pid for root in forest['forest']for group in root['children']for pid in group['partIds']]
            self.assertEqual(len(leaves),len(set(leaves)))
            self.assertEqual(set(leaves),set(p))


if __name__=='__main__':unittest.main()
