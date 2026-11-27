"""Actual packed reciprocal seams must survive a catalogue pose revision."""
import copy,importlib.util,json,tempfile,unittest
from pathlib import Path
from unittest.mock import patch
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location('revision_tissue',ROOT/'scripts/local_catalogue_revision.py')
revision=importlib.util.module_from_spec(spec);spec.loader.exec_module(revision)

class ReciprocalSeams(unittest.TestCase):
    def setUp(self):
        self.tmp=tempfile.TemporaryDirectory();self.addCleanup(self.tmp.cleanup)
        self.root=Path(self.tmp.name);folder=self.root/'public';folder.mkdir()
        payload=bytearray()
        def array(values,dtype):
            a=np.array(values,dtype=dtype);offset=len(payload);payload.extend(a.tobytes());return offset
        positions=array([[0,0,0],[.01,0,0],[0,.01,0]],'<f4')
        normals=array([[0,0,32767]]*3,'<i2')
        signs_a=array([1]*3,'i1');signs_b=array([-1]*3,'i1')
        faces_a=array([0,1,2],'<u4');faces_b=array([0,2,1],'<u4')
        uvs=array([[0,0],[1,0],[0,1]],'<f4')
        self.wrong_normals=array([[0,32767,0]]*3,'<i2')
        (folder/'seam.bin').write_bytes(payload)
        def part(pid,peer,indices,signs):
            return dict(id=pid,chunk=0,positions=positions,vertexCount=3,indices=indices,
                        indexCount=3,normals=normals,normalSigns=signs,uvs=uvs,
                        tissueDrawGroups=[dict(start=0,count=3,neighbor=peer)])
        self.before=dict(chunks=[dict(url='/seam.bin',bytes=len(payload))],
                         parts=[part('heart','cuff',faces_a,signs_a),part('cuff','heart',faces_b,signs_b)])
        self.path=self.root/'atlas.json'

    def guard(self,after,before=None):
        return revision.reciprocal_tissue_seam_guard(before or self.before,after,
                                                    catalogue_path=self.path,root=self.root)

    @staticmethod
    def pose(part):
        part['geometryInstance']=dict(rotation=[.3826834323650898,0,0,.9238795325112867],scale=1.1)
        part['geometryOrigin']=[.1,.2,.3]

    def test_one_sided_heart_pose_rejected_before_install_writes(self):
        after=copy.deepcopy(self.before);self.pose(after['parts'][0])
        with self.assertRaises(revision.TissueSeamRegressionError)as caught:self.guard(after)
        row=caught.exception.evidence['regressions'][0]
        self.assertTrue(row['before']['passed']);self.assertFalse(row['after']['literal_opposed_triangles'])
        cp=self.root/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json'
        cp.parent.mkdir(parents=True);cp.write_bytes(revision.ENC(self.before));original=cp.read_bytes()
        reports=self.root/'reports';reports.mkdir();(reports/'alimentary-visit-paths.json').write_text('[]')
        def mutate(sex,atlas,before):self.pose(atlas['parts'][0]);return dict(changed='heart pose')
        with patch.object(revision,'ROOT',self.root):
            with self.assertRaises(revision.TissueSeamRegressionError):revision.install_revision('seam-fixture',mutate)
        self.assertEqual(cp.read_bytes(),original)
        self.assertFalse((self.root/'.local-models/ta98-review/seam-fixture').exists())

    def test_matched_nontrivial_both_owner_pose_passes(self):
        after=copy.deepcopy(self.before)
        for p in after['parts']:self.pose(p)
        evidence=self.guard(after)
        self.assertTrue(evidence['passed']);self.assertTrue(evidence['allCheckedSeamsPassing'])
        self.assertEqual(evidence['affectedPairs'],1)
        self.assertTrue(evidence['checked'][0]['after']['literal_opposed_normal_words'])

    def test_normal_binding_only_regression_rejected(self):
        after=copy.deepcopy(self.before);after['parts'][1]['normals']=self.wrong_normals
        with self.assertRaises(revision.TissueSeamRegressionError)as caught:self.guard(after)
        actual=caught.exception.evidence['regressions'][0]['after']
        self.assertTrue(actual['literal_opposed_triangles']);self.assertFalse(actual['literal_opposed_normal_words'])

    def test_preexisting_bad_seam_recorded_without_false_pass_and_can_restore(self):
        before=copy.deepcopy(self.before);before['parts'][1]['geometryOrigin']=[0,0,.001]
        after=copy.deepcopy(before);after['parts'][0]['geometryScale']=1.0
        evidence=self.guard(after,before)
        self.assertTrue(evidence['passed']);self.assertFalse(evidence['allCheckedSeamsPassing'])
        self.assertEqual(len(evidence['preExistingBad']),1)
        self.assertFalse(evidence['preExistingBad'][0]['before']['passed'])
        after=copy.deepcopy(before);after['parts'][1]['geometryOrigin']=[0,0,0]
        restored=self.guard(after,before)
        self.assertEqual(len(restored['restored']),1);self.assertTrue(restored['allCheckedSeamsPassing'])

    def test_suppression_visibility_and_names_do_not_decode_or_change_seams(self):
        after=copy.deepcopy(self.before)
        after['parts'][0].update(suppressed=True,defaultHidden=True,name='Renamed heart')
        (self.root/'public/seam.bin').unlink()
        evidence=self.guard(after)
        self.assertTrue(evidence['passed']);self.assertEqual(evidence['affectedPairs'],0)

    def test_removed_or_one_sided_declared_seam_rejected(self):
        for removed in ((0,),(0,1)):
            after=copy.deepcopy(self.before)
            for i in removed:after['parts'][i]['tissueDrawGroups']=[]
            with self.assertRaises(revision.TissueSeamRegressionError):self.guard(after)

    def test_changed_chunk_record_decodes_replacement_normals(self):
        changed=bytearray((self.root/'public/seam.bin').read_bytes())
        # A changed chunk URL must recheck its actual declarations. Equal bytes
        # pass; rebinding only the cuff's normal chart must fail.
        (self.root/'public/replacement.bin').write_bytes(changed)
        after=copy.deepcopy(self.before);after['chunks'][0]['url']='/replacement.bin'
        self.assertTrue(self.guard(after)['allCheckedSeamsPassing'])
        after['parts'][1]['normals']=self.wrong_normals
        with self.assertRaises(revision.TissueSeamRegressionError):self.guard(after)

if __name__=='__main__':unittest.main()
