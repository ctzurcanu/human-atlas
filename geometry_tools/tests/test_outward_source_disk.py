import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.outward_source_disk import select_outward_source_disk

def fixture(center_sign=1,base_sign=-1):
    v=np.array([[0.,0.,0.],[.01,0.,0.],[0.,.01,0.],[0.,0.,.01],[.002,.002,.001]])
    f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3],[0,4,1],[1,4,2],[2,4,0]])
    n=np.tile([0,0,base_sign*32767],(5,1)).astype('<i2');n[4]=[0,0,center_sign*32767]
    return Mesh(v,f,'test:metres',dict(normals=n,normalSigns=np.ones(5,np.int8)))

class OutwardSourceDiskTests(unittest.TestCase):
    def test_redundant_inward_patch_removed_without_changing_the_outer_tetrahedron(self):
        m=fixture();keep,field,index=select_outward_source_disk(m,[[0],[4,5,6]],representation='outer-source-envelope',minimum_outward_cosine=.1)
        np.testing.assert_array_equal(keep,[0,1,2,3]);self.assertEqual(field['removedAlternative'],1);self.assertTrue(index['sourceOpenBoundaryEdgesLiteral'])
    def test_two_consistent_choices_remain_ambiguous(self):
        with self.assertRaisesRegex(ValueError,'unique'):select_outward_source_disk(fixture(-1),[[0],[4,5,6]],representation='outer-source-envelope',minimum_outward_cosine=.1)
    def test_two_bad_choices_are_not_accepted(self):
        with self.assertRaisesRegex(ValueError,'unique'):select_outward_source_disk(fixture(1,1),[[0],[4,5,6]],representation='outer-source-envelope',minimum_outward_cosine=.1)
    def test_a_layer_or_leaflet_cannot_be_classified_by_this_rule(self):
        with self.assertRaisesRegex(ValueError,'classification'):select_outward_source_disk(fixture(),[[0],[4,5,6]],representation='valve-leaflet',minimum_outward_cosine=.1)

if __name__=='__main__':unittest.main()
