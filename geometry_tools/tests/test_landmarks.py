from pathlib import Path
import copy
import hashlib
import tempfile
import unittest
import numpy as np

from geometry_tools.model import AtlasReader
from geometry_tools.landmarks import validate_catalogue


class PointReferenceTests(unittest.TestCase):
    def fixture(self,folder):
        root=Path(folder);path=root/'.local-models/host.bin';path.parent.mkdir()
        path.write_bytes(np.array([[0,0,0],[1,0,0],[0,1,0]],'<f4').tobytes()+np.array([0,1,2],'<u4').tobytes())
        part=dict(id='wall',chunk=0,positions=0,indices=36,vertexCount=3,indexCount=3,tissueDrawGroups=[dict(start=0,count=3,neighbor='ligament')])
        ref=dict(representation='point-landmark',partId='wall',point=[.25,.25,0],face=0,barycentric=[.5,.25,.25],sourceBankSha256=hashlib.sha256(path.read_bytes()).hexdigest(),sourcePlanSha256='a'*64)
        atlas=dict(sex='male',parts=[part],chunks=[dict(url='/local-models/host.bin')],concepts=[dict(id='apex',elements=['wall'],ta98Term='A08.3.01.002')],
                   finiteUrinaryIntegration=dict(apicalAttachment=dict(wallPartId='wall',partId='ligament',sourcePlanSha256='a'*64),termGroups=[dict(term='A08.3.01.002',conceptId='apex',partIds=['wall'],pointReference=ref)]))
        return atlas,AtlasReader(root),path

    def test_reference_is_bound_to_encoded_geometry_without_an_extra_solid(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder);before=copy.deepcopy(atlas['parts'])
            checks=validate_catalogue(atlas,reader=reader,frame='test')
            self.assertEqual(checks[0]['roundtrip_error_m'],0.)
            self.assertFalse(checks[0]['new_tissue_solid']);self.assertEqual(atlas['parts'],before)

    def test_host_translation_cannot_leave_a_stale_point_behind(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder);atlas['parts'][0]['geometryOrigin']=[0,0,.01]
            with self.assertRaisesRegex(ValueError,'encoded host facet'):validate_catalogue(atlas,reader=reader,frame='test')

    def test_changed_source_bank_is_rejected_even_with_valid_bounds(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,path=self.fixture(folder);data=bytearray(path.read_bytes());data[:4]=np.array([.1],'<f4').tobytes();path.write_bytes(data)
            with self.assertRaisesRegex(ValueError,'source bank changed'):validate_catalogue(atlas,reader=reader,frame='test')

    def test_point_cannot_change_to_a_different_support_interface(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder);atlas['parts'][0]['tissueDrawGroups'][0]['neighbor']='another organ'
            with self.assertRaisesRegex(ValueError,'wall-to-ligament'):validate_catalogue(atlas,reader=reader,frame='test')


class DiaphragmReferenceTests(unittest.TestCase):
    def fixture(self,folder):
        atlas,reader,path=PointReferenceTests().fixture(folder)
        atlas.pop('finiteUrinaryIntegration')
        terms=[f'A04.4.02.{i:03d}'for i in range(1,18)]
        refs={10,11,14,15,16}
        atlas['parts']=[dict(atlas['parts'][0],id=term,ta98Term=term)for i,term in enumerate(terms,1)if i!=1 and i not in refs]
        ids=[p['id']for p in atlas['parts']]
        bank=hashlib.sha256(path.read_bytes()).hexdigest()
        def boundary(name):
            return dict(name=name,closed=True,points=[[.25,0,0],[.25,.25,0],[0,.25,0]],
                samplePartIds=[terms[16]]*3,faceIds=[0]*3,barycentric=[[.75,.25,0],[.5,.25,.25],[.75,0,.25]],sourceBankSha256=bank)
        groups=[]
        for i,term in enumerate(terms,1):
            group=dict(term=term,conceptId='concept:'+term,partIds=ids if i==1 else[terms[16]]if i in refs else[term],
                       representation='whole-organ'if i==1 else'boundary-reference'if i in refs else'finite-tissue')
            if i in refs:group['boundaries']=[boundary('left'),boundary('right')]if i in[15,16]else[boundary('rim')]
            if i==2:
                group['representation']='tissue-assembly'
                group['partIds']=[pid for pid in ids if pid in terms[1:7]]
            groups.append(group)
        atlas['concepts']=[dict(id=g['conceptId'],elements=g['partIds'],ta98Term=g['term'],ta98Kind='exact')for g in groups]
        atlas['finiteDiaphragmIntegration']=dict(schema='diaphragm-complete-local-v7',partIds=ids,bankSha256=bank,sourcePlanSha256='a'*64,termGroups=groups)
        return atlas,reader,path

    def test_all_openings_and_paired_triangles_replay_on_actual_encoded_tissue(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder);before=copy.deepcopy(atlas['parts'])
            checks=validate_catalogue(atlas,reader=reader,frame='test')
            self.assertEqual(len(checks),7)
            self.assertTrue(all(c['passed']and not c['new_tissue_solid']for c in checks))
            self.assertTrue(all(Path(c['bank_path']).is_file()for c in checks))
            self.assertEqual(max(c['maximum_roundtrip_error_m']for c in checks),0)
            self.assertEqual(atlas['parts'],before)

    def test_missing_terms_and_wrong_physical_materials_are_rejected(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder)
            atlas['finiteDiaphragmIntegration']['termGroups'].pop()
            with self.assertRaisesRegex(ValueError,'seventeen'):validate_catalogue(atlas,reader=reader,frame='test')
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder)
            next(p for p in atlas['parts'] if p['ta98Term']=='A04.4.02.008')['ta98Term']='A04.4.02.017'
            with self.assertRaisesRegex(ValueError,'own anatomical term'):validate_catalogue(atlas,reader=reader,frame='test')

    def test_camera_independent_source_translation_invalidates_a_stale_rim(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder)
            atlas['parts'][-1]['geometryOrigin']=[0,0,.001]
            with self.assertRaisesRegex(ValueError,'encoded host facet'):validate_catalogue(atlas,reader=reader,frame='test')

    def test_source_bank_and_outside_barycentrics_are_rejected(self):
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,path=self.fixture(folder)
            data=bytearray(path.read_bytes());data[:4]=np.array([.1],'<f4').tobytes();path.write_bytes(data)
            with self.assertRaisesRegex(ValueError,'reference bank is stale'):validate_catalogue(atlas,reader=reader,frame='test')
        with tempfile.TemporaryDirectory()as folder:
            atlas,reader,_=self.fixture(folder)
            atlas['finiteDiaphragmIntegration']['termGroups'][10]['boundaries'][0]['barycentric'][0]=[-.1,.5,.6]
            with self.assertRaisesRegex(ValueError,'finite barycentrics'):validate_catalogue(atlas,reader=reader,frame='test')


if __name__=='__main__':unittest.main()
