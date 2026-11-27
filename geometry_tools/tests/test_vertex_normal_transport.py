import unittest
import numpy as np
from geometry_tools.algorithms.surface_normals import transport_vertex_words

class VertexNormalTransportTest(unittest.TestCase):
    def test_affine_shear_and_literal_disconnected_star(self):
        p=np.array([[0.,0,0],[1.,0,0],[0,1.,0],[3.,0,0],[4.,0,0],[3.,1,0]])
        f=np.array([[0,1,2],[3,4,5]]);w=np.tile([0,0,32767],(6,1)).astype('<i2')
        q=p.copy();q[:3,2]=.2*p[:3,0]
        result,proof=transport_vertex_words(p,q,f,w,maximum_facet_rotation_degrees=20)
        expected=np.array([-.2,0,1.]);expected/=np.linalg.norm(expected)
        np.testing.assert_allclose(result[:3]/32767,np.tile(expected,(3,1)),atol=1/32767)
        np.testing.assert_array_equal(result[3:],w[3:]);self.assertEqual(proof['changed_vertex_fields'],3)

    def test_rotation_bound_rejects(self):
        p=np.array([[0.,0,0],[1.,0,0],[0,1.,0]]);q=p.copy();q[1,2]=1
        with self.assertRaises(ValueError):
            transport_vertex_words(p,q,[[0,1,2]],np.tile([0,0,32767],(3,1)),maximum_facet_rotation_degrees=20)
