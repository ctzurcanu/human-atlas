import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.projected_clearance import directional_clearance


class DirectionalClearanceTests(unittest.TestCase):
    def fixture(self):
        a=Mesh(np.array([[0.,0.,.01],[1.,0.,.02],[0.,1.,.01]]),np.array([[0,1,2]]),'test-metres')
        b=Mesh(np.array([[.2,.2,0.],[.3,.2,0.],[.2,.3,0.]]),np.array([[0,1,2]]),'test-metres')
        return a,b

    def test_extremum_in_projection_interior_not_at_a_moving_vertex(self):
        a,b=self.fixture();r=directional_clearance(a,b)
        self.assertAlmostEqual(r['maximum_forward_excess_m'],.013)
        self.assertTrue(np.allclose(r['moving_witness_m'],[.3,.2,.013]))
        shifted=Mesh(a.vertices-[0,0,.015],a.faces,a.frame)
        self.assertAlmostEqual(directional_clearance(shifted,b)['maximum_forward_excess_m'],-.002)

    def test_rotated_direction_preserves_the_same_clearance(self):
        a,b=self.fixture();rotation=np.array([[0,0,1],[0,1,0],[-1,0,0]])
        aa=Mesh(a.vertices@rotation.T,a.faces,a.frame);bb=Mesh(b.vertices@rotation.T,b.faces,b.frame)
        self.assertAlmostEqual(directional_clearance(aa,bb,direction=(1.,0.,0.))['maximum_forward_excess_m'],.013)

    def test_disjoint_projected_triangles_are_not_called_a_join(self):
        a,b=self.fixture();b=Mesh(b.vertices+[2,0,0],b.faces,b.frame)
        self.assertFalse(directional_clearance(a,b)['projection_overlaps'])


if __name__=='__main__':unittest.main()
