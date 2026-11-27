import unittest
import numpy as np
from geometry_tools.algorithms.compact_translation import CompactTranslation,bound_composition


class CompactTranslationTests(unittest.TestCase):
    def setUp(self):
        self.field = CompactTranslation(np.array([.03,1.2,-.02]),
                                        np.array([.03,.07,.025]),
                                        np.array([-.00035,0.,.0001]),
                                        'model:test:assembled-metres')
        self.rng = np.random.default_rng(194)

    def test_analytic_derivative_matches_independent_finite_difference(self):
        points = self.field.center+self.rng.uniform(-.6,.6,(120,3))*self.field.radii
        numeric = np.empty((len(points),3,3))
        for axis in range(3):
            delta = np.eye(3)[axis]*1e-7
            numeric[:,:,axis] = (self.field.map(points+delta)-self.field.map(points-delta))/(2e-7)
        np.testing.assert_allclose(self.field.jacobians(points),numeric,atol=2e-9,rtol=1e-7)

    def test_global_bound_and_inverse_are_consistent(self):
        points = self.field.center+self.rng.uniform(-1.2,1.2,(10000,3))*self.field.radii
        J = self.field.jacobians(points)
        smallest = np.linalg.svd(J,compute_uv=False)[:,-1]
        self.assertGreaterEqual(smallest.min(),self.field.minimum_stretch_bound-1e-12)
        self.assertGreaterEqual(np.linalg.det(J).min(),self.field.minimum_stretch_bound-1e-12)
        inverse = self.field.inverse_transposes(points).transpose(0,2,1)
        np.testing.assert_allclose(inverse@J,np.broadcast_to(np.eye(3),J.shape),atol=5e-16)

    def test_transported_normals_remain_perpendicular_to_transported_tangents(self):
        points = self.field.center+self.rng.uniform(-.6,.6,(100,3))*self.field.radii
        first,second = self.rng.normal(size=(2,len(points),3))
        normal = np.cross(first,second)
        J = self.field.jacobians(points)
        n = self.field.transform_normals(points,normal)
        for tangent in (first,second):
            moved = np.einsum('nij,nj->ni',J,tangent)
            np.testing.assert_allclose(np.einsum('ni,ni->n',n,moved),0.,atol=2e-14)
        np.testing.assert_allclose(self.field.transform_normals(points,-normal),-n,atol=1e-15)

    def test_compact_boundary_and_source_preservation(self):
        points = np.array([self.field.center,
                           self.field.center+self.field.radii*np.array([1.,0.,0.]),
                           self.field.center+self.field.radii*3])
        before = points.copy();mapped = self.field.map(points)
        np.testing.assert_array_equal(points,before)
        np.testing.assert_allclose(mapped[0],points[0]+self.field.displacement,atol=0)
        np.testing.assert_array_equal(mapped[1:],points[1:])
        np.testing.assert_array_equal(self.field.jacobians(points[1:]),np.broadcast_to(np.eye(3),(2,3,3)))

    def test_noninjective_parameters_and_invalid_normals_are_rejected(self):
        with self.assertRaises(ValueError):
            CompactTranslation(np.zeros(3),np.ones(3)*.01,np.array([.1,0,0]),'test')
        with self.assertRaises(ValueError):
            self.field.transform_normals(np.zeros((1,3)),np.zeros((1,3)))

    def test_complete_box_bounds_cover_positions_and_local_derivatives(self):
        low=self.field.center+np.array([[-.01,-.02,-.01],[.02,.01,.01],[.08,.08,.08]])
        high=low+np.array([.012,.013,.012])
        moved_low,moved_high,stretch=self.field.bound_boxes(low,high)
        for i in range(len(low)):
            points=self.rng.uniform(low[i],high[i],(2000,3))
            mapped=self.field.map(points)
            self.assertTrue(np.all(mapped>=moved_low[i])and np.all(mapped<=moved_high[i]))
            singular=np.linalg.svd(self.field.jacobians(points),compute_uv=False)[:,-1]
            self.assertGreaterEqual(singular.min(),stretch[i]-1e-12)
        self.assertEqual(stretch[-1],1.)

    def test_composed_box_certificate_is_not_a_corner_sampling_test(self):
        second=CompactTranslation(self.field.center+np.array([.08,0,0]),
                                  self.field.radii,self.field.displacement,self.field.frame)
        low=np.array([self.field.center-self.field.radii,second.center-second.radii])
        high=np.array([self.field.center+self.field.radii,second.center+second.radii])
        _,_,bound=bound_composition([self.field,second],low,high)
        for i in range(2):
            points=self.rng.uniform(low[i],high[i],(2000,3))
            J=second.jacobians(self.field.map(points))@self.field.jacobians(points)
            self.assertGreaterEqual(np.linalg.svd(J,compute_uv=False)[:,-1].min(),bound[i]-1e-12)
        self.assertGreater(bound.min(),self.field.minimum_stretch_bound*second.minimum_stretch_bound)

    def test_directional_bounds_distinguish_expansion_and_shear_from_contraction(self):
        for center,displacement in (([-.25,0,0],[.1,0,0]),([.25,0,0],[0,.1,0])):
            field=CompactTranslation(np.zeros(3),np.ones(3),np.array(displacement),'test')
            low=np.array([center])-.001;high=np.array([center])+.001
            _,_,bound=field.bound_boxes(low,high)
            points=self.rng.uniform(low[0],high[0],(5000,3))
            actual=np.linalg.svd(field.jacobians(points),compute_uv=False)[:,-1]
            self.assertGreaterEqual(actual.min(),bound[0]-1e-12)
            self.assertGreater(bound[0],field.minimum_stretch_bound+.08)


if __name__ == '__main__':
    unittest.main()
