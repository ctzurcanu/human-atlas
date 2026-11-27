import importlib.util,json,unittest
from pathlib import Path

BASE=Path('.local-models/geometry-workbench/biodigital-male-circulation-v1')
spec=importlib.util.spec_from_file_location('explicit_labial_host','scripts/verify-labial-host-course.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class LabialHostCourseTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.c=json.loads((BASE/'labial-intramuscular-course-v1860/contract.json').read_text())
        cls.r=json.loads((BASE/'labial-intramuscular-course-v1860/report.json').read_text())
        cls.g=json.loads((BASE/'labial-intramuscular-native-v1861/report.json').read_text())
        cls.context=json.loads((BASE/'head-fixed-context-v1665/report.json').read_text())['rows']

    def test_actual_full_sources_preserve_calibre_and_all_end_openings(self):
        self.assertTrue(module.verify_candidate(self.c,self.r)['complete'])

    def test_documented_host_contact_is_scoped_to_actual_lip_geometry(self):
        q=module.qualify_host_events(self.c,self.r['rows'],self.g['newNativeOrOtherBloodPairs'],self.context)
        self.assertEqual(len(q),166);self.assertTrue(all(x['qualified']for x in q))

    def test_bone_and_neighbouring_blood_never_gain_host_exemptions(self):
        event=self.g['newNativeOrOtherBloodPairs'][0]
        for owner in ['ZA:Upper lateral incisor.r:Teeth','BDM:male:human_20_male_cardiovascular_system-right_mental_artery_ID']:
            with self.subTest(owner=owner),self.assertRaises(AssertionError):
                module.qualify_host_events(self.c,self.r['rows'],[dict(event,contextOwner=owner)],self.context)

    def test_opposite_side_muscle_is_not_a_host_for_the_changed_lip_segment(self):
        event=self.g['newNativeOrOtherBloodPairs'][0];other=next(x for x in module.HOSTS if x!=event['contextOwner'])
        with self.assertRaises(AssertionError):module.qualify_host_events(self.c,self.r['rows'],[dict(event,contextOwner=other)],self.context)
