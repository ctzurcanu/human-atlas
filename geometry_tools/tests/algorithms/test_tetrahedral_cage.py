import unittest
import numpy as np
from geometry_tools.algorithms.tetrahedral_cage import RegularTetrahedralCage
from geometry_tools.algorithms.selective_plane_refinement import refine_selected_plane


class TetrahedralCageTests(unittest.TestCase):
    def setUp(self):self.c=RegularTetrahedralCage([0,0,0],.01,[4,4,4])

    def test_interpolation_reproduces_source_points_and_affine_fields(self):
        p=np.random.default_rng(7).uniform(.0001,.0299,(100,3));r=self.c.locate(p)
        rebuilt=np.einsum('nk,nki->ni',r['weights'],self.c.nodes[r['node_ids']])
        np.testing.assert_allclose(rebuilt,p,atol=1e-17)
        matrix=np.array([[1.,.1,.2],[.2,.9,0.],[0.,.1,1.2]])
        target=self.c.nodes@matrix.T+[.1,.2,.3]
        got=np.einsum('nk,nki->ni',r['weights'],target[r['node_ids']])
        np.testing.assert_allclose(got,p@matrix.T+[.1,.2,.3],atol=1e-16)

    def test_fixed_boundary_map_is_identity_outside_and_has_positive_elements(self):
        d=np.zeros_like(self.c.nodes);d[~self.c.boundary]=[.0005,-.0003,.0004]
        proof=self.c.check(d);self.assertTrue(proof['passed'])
        p=np.array([[-.001,.01,.01],[.031,.02,.02],[0.,.02,.01]])
        np.testing.assert_array_equal(self.c.map(p,d),p)
        volume=np.abs(np.linalg.det(self.c.source_edges))*np.linalg.det(self.c.jacobians(d))/6
        self.assertAlmostEqual(float(volume.sum()),.03**3,places=18)

    def test_point_jacobians_match_finite_differences_inside_one_tetrahedron(self):
        d=np.zeros_like(self.c.nodes);d[~self.c.boundary]=np.random.default_rng(9).normal(size=((~self.c.boundary).sum(),3))*.0001
        p=np.array([[.015,.012,.017]]);j=self.c.point_jacobians(p,d);finite=np.zeros_like(j)
        for k in range(3):
            eps=np.eye(3)[k]*1e-8
            finite[:,:,k]=(self.c.map(p+eps,d)-self.c.map(p-eps,d))/(2e-8)
        np.testing.assert_allclose(j,finite,atol=3e-10)

    def test_boundary_motion_and_inverted_elements_are_rejected(self):
        d=np.zeros_like(self.c.nodes);d[0]=[.001,0,0]
        self.assertFalse(self.c.check(d)['passed'])
        with self.assertRaises(ValueError):self.c.map([[.015,.012,.017]],d)
        d=np.zeros_like(self.c.nodes);d[self.c.node_ids([1,1,1])]=[.05,0,0]
        self.assertFalse(self.c.check(d)['passed'])

    def test_complete_selective_conformance_makes_surface_mapping_piecewise_affine(self):
        p=np.array([[-.005,.011,.015],[.025,.005,.018],[.020,.028,.003]])
        f=np.array([[0,1,2]])
        for stage in ['box','grid','diagxy','diagxz','diagyz']:
            for normal,coordinate,scope in self.c.conformance_planes(stage):
                active=self.c.active_faces(p,f,scope)
                r=refine_selected_plane(p,f,normal=normal,coordinate=coordinate,origin=self.c.minimum,active_faces=active)
                p,f=r['positions'],r['faces']
        self.assertTrue(self.c.triangle_conformance(p,f)['passed'])
        d=np.zeros_like(self.c.nodes);d[~self.c.boundary]=[.0005,-.0003,.0004]
        mapped=self.c.map(p,d);centers=p[f].mean(1)
        np.testing.assert_allclose(self.c.map(centers,d),mapped[f].mean(1),atol=2e-17)

if __name__=='__main__':unittest.main()
