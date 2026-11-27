import copy,importlib.util,json,unittest
from pathlib import Path
from geometry_tools.course_io import sha

BASE=Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')
spec=importlib.util.spec_from_file_location('labial_identity_revision','scripts/prepare-current-labial-mouth.py');module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class LabialIdentityRevisionTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.source=BASE/'blood-installed-v1865/circulation-atlas.json'
        cls.old=json.loads(cls.source.read_text());cls.new=json.loads((BASE/'blood-installed-v1876/circulation-atlas.json').read_text())

    def check(self,value):
        return module.verify_identity_revision(self.old,value,sha(self.source),'blood-installed-v1876')

    def test_actual_names_follow_fixed_native_and_facial_artery_hemispheres(self):
        self.assertTrue(self.check(self.new))

    def test_the_naming_revision_cannot_hide_a_skeleton_change(self):
        changed=copy.deepcopy(self.new);bone=next(p for p in changed['parts']if p['system']=='skeletal');bone['bounds'][0][0]+=.001
        with self.assertRaises(AssertionError):self.check(changed)

    def test_no_other_concept_can_be_renamed(self):
        changed=copy.deepcopy(self.new);other=next(c for c in changed['concepts']if c['id']not in {x['partId']for x in module.CORRECTIONS});other['name']='Unexpected replacement'
        with self.assertRaises(AssertionError):self.check(changed)

    def test_corrected_names_cannot_be_swapped_back_under_the_same_review(self):
        changed=copy.deepcopy(self.new);p=next(p for p in changed['parts']if p['id']==module.CORRECTIONS[0]['partId']);p['name']=module.CORRECTIONS[0]['before']
        with self.assertRaises(AssertionError):self.check(changed)
