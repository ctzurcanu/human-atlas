"""Regress stale suppressed-parent reuse after a source-region installation."""
import tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.model import AtlasReader
from geometry_tools.current_source_carrier import reconstruct_current_carrier


class CurrentSourceCarrier(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        folder=self.root/'.local-models';folder.mkdir()
        self.vertices=np.array([[0,0,0],[.01,0,0],[.01,.01,0],[0,.01,0]],'<f4')
        self.current=self.vertices+np.array([.02,.03,0],'<f4')
        normals=np.tile([0,0,32767],(4,1)).astype('<i2');signs=np.ones(4,'i1')
        faces=np.array([[0,1,2],[0,2,3]],'<u4')
        for name,points in [('old',self.vertices),('current',self.current)]:
            (folder/(name+'.bin')).write_bytes(points.tobytes()+normals.tobytes()+signs.tobytes()+faces.tobytes())
        base=dict(vertexCount=4,positions=0,normals=48,normalSigns=72,indices=76,indexCount=6)
        self.atlas=dict(chunks=[dict(url='/local-models/old.bin'),dict(url='/local-models/current.bin')],defaultHidden=['parent'],
                        parts=[dict(base,id='parent',chunk=0,suppressed=True),
                               dict(base,id='first',chunk=1,indexCount=3,sourceDerivedRegion=dict(sourceFaceIds=[0])),
                               dict(base,id='second',chunk=1,indices=88,indexCount=3,sourceDerivedRegion=dict(sourceFaceIds=[1]))])
        self.reader=AtlasReader(self.root);self.source=self.reader.read_part(self.atlas,self.atlas['parts'][0],frame='test:metres')
        self.row=dict(activePartIds=['first','second'])

    def tearDown(self):self.reader.maps.clear();self.temp.cleanup()

    def test_active_partitions_replace_stale_parent_geometry(self):
        actual,projections=reconstruct_current_carrier(self.atlas,self.source,self.row,self.reader,catalogue_path=None)
        np.testing.assert_array_equal(actual.vertices,self.current.astype(float))
        np.testing.assert_array_equal(actual.faces,self.source.faces)
        self.assertEqual(len(projections),2)

    def test_missing_partition_cannot_claim_a_complete_source(self):
        with self.assertRaises(AssertionError):
            reconstruct_current_carrier(self.atlas,self.source,dict(activePartIds=['first']),self.reader,catalogue_path=None)

    def test_suppressed_source_is_never_an_active_partition(self):
        with self.assertRaises(AssertionError):
            reconstruct_current_carrier(self.atlas,self.source,dict(activePartIds=['parent']),self.reader,catalogue_path=None)

if __name__=='__main__':unittest.main()
