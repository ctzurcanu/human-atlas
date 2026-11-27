import copy,importlib.util,unittest
from pathlib import Path

spec=importlib.util.spec_from_file_location('local_revision',Path(__file__).resolve().parents[1]/'scripts/local_catalogue_revision.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

class VisitPreservation(unittest.TestCase):
    def setUp(self):
        self.atlas=dict(parts=[dict(id='tube',chunk=0,positions=0,indices=96,geometryOrigin=[0,1,0]),dict(id='cord',chunk=1)],chunks=[dict(url='/tube.bin',bytes=120),dict(url='/cord.bin',bytes=120)])
        self.visit=dict(bindings=[dict(id='tube')],sourceCatalogueSha256='stale')
    def test_unrelated_organ_keeps_literal_old_route(self):
        after=copy.deepcopy(self.atlas);after['parts'][1]['positions']=24
        self.assertTrue(module.stale_visit_geometry_unchanged(self.atlas,after,self.visit))
    def test_route_part_geometry_pose_or_visibility_change_rejects(self):
        for key,value in [('positions',12),('geometryOrigin',[0,1.01,0]),('suppressed',True)]:
            after=copy.deepcopy(self.atlas);after['parts'][0][key]=value
            self.assertFalse(module.stale_visit_geometry_unchanged(self.atlas,after,self.visit))
    def test_replaced_bank_or_missing_part_rejects(self):
        after=copy.deepcopy(self.atlas);after['chunks'][0]['url']='/replacement.bin'
        self.assertFalse(module.stale_visit_geometry_unchanged(self.atlas,after,self.visit))

    def tongue_retirement(self):
        before=copy.deepcopy(self.atlas);before['parts'][0]['conceptId']='native'
        after=copy.deepcopy(before);after['parts'][0]['suppressed']=True
        after['parts'].append(dict(id='new-tongue',chunk=0,indexCount=12))
        after['concepts']=[dict(id='native',elements=['new-tongue']),dict(id='whole',elements=['new-tongue'])]
        after['finiteTongueIntegration']=dict(schema='tongue-local/v1',legacyWholeAlias=True,sourcePartId='tube',partIds=['new-tongue'],wholeConceptId='whole')
        return before,after

    def test_literal_retired_tongue_keeps_stale_route_without_certification(self):
        before,after=self.tongue_retirement();original=copy.deepcopy(self.visit)
        self.assertTrue(module.stale_visit_geometry_unchanged(before,after,self.visit))
        self.assertEqual(self.visit,original)
        self.assertEqual(self.visit['sourceCatalogueSha256'],'stale')

    def test_retirement_cannot_hide_geometry_or_alias_change(self):
        for change in ('pose','alias','hidden-replacement','missing-source'):
            before,after=self.tongue_retirement()
            if change=='pose':after['parts'][0]['geometryOrigin']=[0,2,0]
            elif change=='alias':after['concepts'][0]['elements']=['tube']
            elif change=='hidden-replacement':after['parts'][-1]['suppressed']=True
            else:after['parts'].pop(0)
            self.assertFalse(module.stale_visit_geometry_unchanged(before,after,self.visit),change)
        after=copy.deepcopy(self.atlas);after['parts'].pop(0)
        self.assertFalse(module.stale_visit_geometry_unchanged(self.atlas,after,self.visit))

if __name__=='__main__':unittest.main()
