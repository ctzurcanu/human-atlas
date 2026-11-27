import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.occupancy_classifier import close_literal_rims
from geometry_tools.audit import topology


class OccupancyClassifierTests(unittest.TestCase):
    def test_open_square_tube_retains_sides_and_has_only_analysis_caps(self):
        v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,2],[1,0,2],[1,1,2],[0,1,2]],float)
        f=np.array([[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        source=Mesh(v,f,'test:metres');closed,proof=close_literal_rims(source);self.assertTrue(proof['analysisOnly']);self.assertFalse(proof['deliveryPermitted']);self.assertEqual(len(proof['literalRims']),2);self.assertEqual(proof['classifierFaces'],8);self.assertAlmostEqual(topology(closed)['oriented_triangle_integral_m3'],2.)
        np.testing.assert_array_equal(closed.triangles[:len(f)],source.triangles)

    def test_nonconvex_rim_is_rejected(self):
        p=np.array([[0,0,0],[2,0,0],[.2,.2,0],[0,2,0]],float);v=np.r_[p,p+[0,0,1]];f=[]
        for i in range(4):j=(i+1)%4;f.extend([[i,j,j+4],[i,j+4,i+4]])
        with self.assertRaises(ValueError):close_literal_rims(Mesh(v,np.array(f),'test:metres'))
