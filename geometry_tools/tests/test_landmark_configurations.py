"""Clinical opening witnesses retain the sex-specific TA98 identity."""
import copy
from types import SimpleNamespace
import unittest
import numpy as np
from geometry_tools.landmarks import validate_catalogue


class ConfigurationTests(unittest.TestCase):
    def fixture(self,sex):
        term,opening={'male':('A09.4.02.002','A08.3.01.027M'),
                      'female':('A09.2.03.002','A08.3.01.028F')}[sex]
        sha='a'*64;pid='lining';parent='opening';cid='configuration'
        boundary=dict(name='neck',partIds=[pid],sourcePlanSha256=sha,
            points=[[.25,.25,0]],samplePartIds=[pid],faceIds=[0],barycentric=[[.5,.25,.25]])
        atlas=dict(sex=sex,parts=[dict(id=pid,tissueDrawGroups=[dict(start=0,count=3)])],
            concepts=[dict(id=cid,ta98Term=term,partOf=parent,elements=[pid])],
            finiteUrinaryIntegration=dict(openingPlanSha256=sha,termGroups=[
                dict(term=opening,conceptId=parent,partIds=[pid],openingBoundaries=[boundary]),
                dict(term=term,conceptId=cid,partIds=[pid],openingConfiguration=dict(
                    state='voiding',openingTerm=opening,sourcePlanSha256=sha,sourceBankSha256=sha))]))
        mesh=SimpleNamespace(triangles=np.array([[[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]]),
            faces=np.array([[0,1,2]]),source=dict(bank_path='fixture.bin',bank_sha256=sha))
        class Reader:
            def read_part(self,atlas,part,*,frame):
                self.frame=frame
                return mesh
        return atlas,Reader()

    def test_both_sexes_replay_their_literal_clinical_rim(self):
        for sex,model in [('male','local-ta98'),('female','local-female-ta98')]:
            with self.subTest(sex=sex):
                atlas,reader=self.fixture(sex);result=validate_catalogue(atlas,reader=reader)
                self.assertEqual(len(result),1);self.assertTrue(result[0]['passed'])
                self.assertEqual(result[0]['maximum_roundtrip_error_m'],0)
                self.assertEqual(reader.frame,'model:'+model+':assembled-metres')

    def test_opposite_sex_and_filling_do_not_pass(self):
        for sex in ['male','female']:
            atlas,reader=self.fixture(sex)
            wrong=copy.deepcopy(atlas);wrong['sex']='female'if sex=='male'else'male'
            with self.assertRaises(ValueError):validate_catalogue(wrong,reader=reader)
            atlas['finiteUrinaryIntegration']['termGroups'][-1]['openingConfiguration']['state']='filling'
            with self.assertRaises(ValueError):validate_catalogue(atlas,reader=reader)

    def test_source_facet_replay_is_still_required(self):
        atlas,reader=self.fixture('female')
        atlas['finiteUrinaryIntegration']['termGroups'][0]['openingBoundaries'][0]['points'][0][2]=.001
        with self.assertRaises(ValueError):validate_catalogue(atlas,reader=reader)


if __name__=='__main__':unittest.main()
