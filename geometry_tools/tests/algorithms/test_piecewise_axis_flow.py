import unittest
import numpy as np
from geometry_tools.algorithms.piecewise_axis_flow import conformed_axis_flow,transport_axis_corner_normals
from geometry_tools.algorithms.shared_plane_refinement import refine_at_axis_plane


class PiecewiseAxisFlowTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]])
        self.f=np.array([[0,1,2],[0,2,3]])
        self.spec=dict(axis=0,knots=[-1.,0.,1.],targets=[[-1.,0.,0.],[.2,.1,.3],[.8,.2,-.1]])

    def test_unrefined_crossing_is_rejected(self):
        with self.assertRaises(ValueError):conformed_axis_flow(self.p,self.f,**self.spec)

    def test_shared_refinement_makes_exact_affine_patches(self):
        r=refine_at_axis_plane(self.p,self.f,axis=0,coordinate=0.)
        m=conformed_axis_flow(r['positions'],r['faces'],**self.spec)
        self.assertAlmostEqual(m['evidence']['analytic_determinant_lower_bound'],.6)
        self.assertEqual(len(m['positions']),len(r['positions']))
        face=m['positions'][r['faces']]
        self.assertTrue(np.all(np.linalg.norm(np.cross(face[:,1]-face[:,0],face[:,2]-face[:,0]),axis=1)>0))
        ids=np.flatnonzero(r['positions'][:,0]==0.)
        np.testing.assert_allclose(m['positions'][ids,0],.2)

    def test_inverse_transpose_transport_matches_explicit_matrices(self):
        slopes=np.array([[1.2,.1,.4],[.7,-.3,.2]])
        normals=np.random.default_rng(8).normal(size=(2,3,3))
        got=transport_axis_corner_normals(normals,slopes,axis=0)
        matrix=np.broadcast_to(np.eye(3),(2,3,3)).copy();matrix[:,:,0]=slopes
        expected=np.einsum('tji,tkj->tki',np.linalg.inv(matrix),normals)
        np.testing.assert_allclose(got,expected,atol=3e-16)

    def test_reversing_primary_coordinate_is_rejected(self):
        spec=dict(self.spec);spec['targets']=[[-1.,0.,0.],[.2,.1,.3],[.1,.2,-.1]]
        r=refine_at_axis_plane(self.p,self.f,axis=0,coordinate=0.)
        with self.assertRaises(ValueError):conformed_axis_flow(r['positions'],r['faces'],**spec)

if __name__=='__main__':unittest.main()
