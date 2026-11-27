import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.annular_pipe_envelopes import recover_paired_sheets
from geometry_tools.algorithms.simple_compartment_section import simple_compartment_section,polygon_clearance
from geometry_tools.tests.algorithms.test_root_collar_transport import finite_pipe


class SimpleCompartmentSectionTests(unittest.TestCase):
    def fluid(self):
        p,f,_=finite_pipe(8);labels=np.r_[np.ones(16,int),np.zeros(16,int)]
        r=recover_paired_sheets(p,f,labels)
        return Mesh(r.positions,r.fluid_faces,'fixture:metres')
    def test_an_off_centre_guide_recovers_the_actual_polygon_and_centroid(self):
        m=self.fluid();r=simple_compartment_section(m,[.2,.1,.37],[0,0,1])
        np.testing.assert_allclose(r['centroid_m'],[0,0,.37],atol=1e-14)
        self.assertAlmostEqual(r['area_m2'],4*.8**2*np.sin(np.pi/4),places=13)
        self.assertGreater(r['centroid_clearance_m'],r['query_clearance_m'])
        self.assertEqual(r['caps_created'],0)
    def test_disjoint_components_are_not_joined_by_a_section(self):
        m=self.fluid();p=np.r_[m.vertices,m.vertices+[4,0,0]];f=np.r_[m.faces,m.faces+len(m.vertices)]
        r=simple_compartment_section(Mesh(p,f,m.frame),[4.1,0,.37],[0,0,1])
        np.testing.assert_allclose(r['centroid_m'],[4,0,.37],atol=1e-14)
        self.assertEqual(r['other_disjoint_contours'],1)
    def test_wall_annulus_cannot_be_mistaken_for_a_solid_fluid_section(self):
        p,f,_=finite_pipe(8)
        with self.assertRaisesRegex(ValueError,'Hole-bearing'):
            simple_compartment_section(Mesh(p,f,'fixture:metres'),[.9,0,.37],[0,0,1])
    def test_open_surfaces_and_external_queries_are_explicit(self):
        m=self.fluid()
        with self.assertRaisesRegex(ValueError,'No section'):
            simple_compartment_section(m,[3,0,.37],[0,0,1])
        with self.assertRaisesRegex(ValueError,'open surface'):
            simple_compartment_section(Mesh(m.vertices,m.faces[1:],m.frame),[0,0,.37],[0,0,1])
    def test_literal_repeated_projected_points_do_not_change_the_boundary(self):
        square=np.array([[0.,0.],[1.,0.],[1.,1.],[0.,1.]])
        repeated=square[[0,1,1,2,3,0]]
        self.assertEqual(polygon_clearance(repeated,np.array([.25,.4])),.25)
        self.assertEqual(polygon_clearance(repeated,np.array([1.2,.4])),
                         polygon_clearance(square,np.array([1.2,.4])))
    def test_short_edges_remain_and_a_fully_collapsed_projection_is_rejected(self):
        tiny=np.ldexp(1.,-40)
        square=np.array([[0.,0.],[tiny,0.],[tiny,tiny],[0.,tiny]])
        self.assertEqual(polygon_clearance(square,np.array([tiny/2,tiny/2])),tiny/2)
        with self.assertRaisesRegex(ValueError,'no nonzero projected edges'):
            polygon_clearance(np.zeros((4,2)),np.array([0.,0.]))
    def test_declared_plane_roundoff_bound_preserves_a_near_vertex_section(self):
        # A literal almost-on-plane vertex makes a crossing fraction round to
        # one if the declared classification bound is not applied to the field.
        p=np.array([[0.,1.,1.],[0.,0.,-1e-30],[1.,0.,-1.],[0.,0.,-1.]])
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        m=Mesh(p.copy(),f.copy(),'fixture:metres')
        with self.assertRaisesRegex(ValueError,'not representable'):
            simple_compartment_section(m,[.1,.2,0.],[0.,0.,1.],classification_tolerance_m=0.)
        r=simple_compartment_section(m,[.1,.2,0.],[0.,0.,1.])
        self.assertAlmostEqual(r['area_m2'],.125,places=15)
        self.assertEqual(r['controls_classified_on_plane'],1)
        self.assertFalse(r['source_geometry_changed'])
        np.testing.assert_array_equal(m.vertices,p);np.testing.assert_array_equal(m.faces,f)
    def test_invalid_plane_classification_bounds_are_rejected(self):
        for bound in [-1.,np.inf,np.nan]:
            with self.assertRaisesRegex(ValueError,'nonnegative plane-classification'):
                simple_compartment_section(self.fluid(),[.1,.2,.4],[0.,0.,1.],classification_tolerance_m=bound)


if __name__=='__main__':unittest.main()
