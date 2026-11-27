import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.section_pair_separation import triangle_axis_rows, separated_triangle_axes,profile_radial_axes,directed_triangle_exit_axes
from geometry_tools.algorithms.triangle_cross_candidates import TriangleCrossIndex


class TriangleAxisTests(unittest.TestCase):
    def test_directed_exit_certifies_a_real_triangle_interval(self):
        a,b=self.meshes();axis,amount=directed_triangle_exit_axes(a.triangles,b.triangles,[1.,0.,0.],clearance_m=.01)
        moved=a.triangles+np.array([1.,0.,0.])*(amount[0]+1e-8)
        gap=(moved[0]@axis[0]).min()-(b.triangles[0]@axis[0]).max()
        self.assertGreaterEqual(gap,.01)
        self.assertGreater(axis[0,0],0)

    def test_directed_exit_handles_a_grazing_obstacle_plane(self):
        a,_=self.meshes();b=Mesh(np.array([[.2,-10.,-1.],[.2,10.,-1.],[.2,0.,1.]]),np.array([[0,1,2]]),'fixture')
        axis,amount=directed_triangle_exit_axes(a.triangles,b.triangles,[0.,0.,1.],clearance_m=.001)
        self.assertTrue(np.isfinite(amount).all());self.assertLess(amount[0],2.)
        gap=((a.triangles[0]+[0.,0.,amount[0]+1e-8])@axis[0]).min()-(b.triangles[0]@axis[0]).max()
        self.assertGreaterEqual(gap,.001-1e-12)

    def test_profile_axis_uses_radial_centres_and_not_longitudinal_offset(self):
        ring=np.array([[0.,1.,0.],[0.,-.5,.8660254],[0.,-.5,-.8660254]])
        a=Mesh(ring,np.array([[0,1,2]]),'fixture')
        b=Mesh(ring+np.array([4.,-.2,0.]),a.faces,'fixture')
        axis=profile_radial_axes(a,b,np.zeros(3,int),np.zeros(3,int),[[0,0]],preferred_axis=[0.,1.,0.])
        np.testing.assert_allclose(axis,[[0.,1.,0.]],atol=1e-12)

    def test_protected_core_uses_source_graph_without_fabricating_a_profile(self):
        ring=np.array([[0.,1.,0.],[0.,-.5,.8660254],[0.,-.5,-.8660254]])
        vertices=np.r_[np.array([[-1.,0.,0.]]),ring]
        a=Mesh(vertices,np.array([[0,1,2],[1,2,3]]),'fixture')
        b=Mesh(vertices+np.array([3.,-.2,0.]),a.faces,'fixture')
        station=np.array([0,1,1,1])
        axis=profile_radial_axes(a,b,station,station,[[0,0]],preferred_axis=[0.,1.,0.])
        np.testing.assert_allclose(axis,[[0.,1.,0.]],atol=1e-12)

    def test_coincident_profile_faces_cannot_reverse_the_declared_escape_side(self):
        ring=np.array([[0.,1.,0.],[0.,-.5,.8660254],[0.,-.5,-.8660254]])
        a=Mesh(np.r_[ring,ring+[2.,0.,0.]],np.array([[0,1,3],[1,4,3]]),'fixture')
        b=Mesh(np.r_[ring,ring+[2.,0.,0.]],a.faces,'fixture')
        station=np.array([0,0,0,1,1,1])
        axes=profile_radial_axes(a,b,station,station,[[0,1],[1,0]],preferred_axis=[0.,1.,0.])
        self.assertTrue(np.all(axes[:,1]>0))
        np.testing.assert_allclose(axes,[[0.,1.,0.],[0.,1.,0.]],atol=1e-12)

    def meshes(self):
        a = Mesh(np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]]), np.array([[0, 1, 2]]), 'fixture')
        b = Mesh(np.array([[0., 0., -1.], [0., 0., 1.], [0., 1., 0.]]), np.array([[0, 1, 2]]), 'fixture')
        return a, b

    def test_nine_constraints_certify_the_entire_projected_intervals(self):
        a, b = self.meshes()
        wa, wb = np.ones((3, 1)), np.zeros((3, 1))
        w, n, rhs = triangle_axis_rows(a, b, wa, wb, [[0, 0]], [2, 0, 0], gap_m=.1)
        self.assertEqual(w.shape, (9, 1))
        shift = np.array([[.1001, 0., 0.]])
        self.assertTrue(np.all(np.einsum('ni,ni->n', w @ shift, n) >= rhs))
        moved = a.vertices + wa @ shift
        self.assertGreaterEqual(moved[:, 0].min() - b.vertices[:, 0].max(), .1)
        shift[0, 0] = .0999
        self.assertFalse(np.all(np.einsum('ni,ni->n', w @ shift, n) >= rhs))

    def test_relative_motion_cancels_a_common_translation(self):
        a, b = self.meshes()
        weights = np.ones((3, 1))
        w, _, _ = triangle_axis_rows(a, b, weights, weights, [[0, 0]], [0, 0, 1], gap_m=0.)
        np.testing.assert_array_equal(w, 0)

    def test_rejects_invalid_or_unbound_inputs(self):
        a, b = self.meshes()
        for pairs, axis, gap in [([[0, 0]], [0, 0, 0], 0), ([[1, 0]], [1, 0, 0], 0),
                                 ([[0., 0.]], [1, 0, 0], 0), ([[0, 0]], [1, 0, 0], -.1)]:
            with self.subTest(pairs=pairs, axis=axis, gap=gap), self.assertRaises(ValueError):
                triangle_axis_rows(a, b, np.ones((3, 1)), np.zeros((3, 1)), pairs, axis, gap_m=gap)

    def test_separated_coplanar_triangles_keep_a_real_interval_gap(self):
        a,b=self.meshes();ta=a.triangles;tb=ta+np.array([2.,0.,0.])
        axis,gap=separated_triangle_axes(ta,tb)
        self.assertGreater(gap[0],0)
        actual=(ta[0]@axis[0]).min()-(tb[0]@axis[0]).max()
        self.assertAlmostEqual(actual,gap[0])
        _,n,_=triangle_axis_rows(a,b,np.ones((3,1)),np.zeros((3,1)),[[0,0]],axis,gap_m=0)
        np.testing.assert_allclose(n,np.broadcast_to(axis,(9,3)))

    def test_coordinate_preservation_axis_is_valid_for_the_actual_triangles(self):
        a,_=self.meshes();first=np.repeat(a.triangles,2,axis=0);second=first+np.array([[2.,0.,0.],[0.,0.,-2.]])[:,None]
        axes,gaps=separated_triangle_axes(first,second,prefer_coordinates=True)
        self.assertTrue(np.all(gaps>0))
        for ta,tb,axis,gap in zip(first,second,axes,gaps):
            self.assertAlmostEqual((ta@axis).min()-(tb@axis).max(),gap)

    def test_padded_boxes_include_possible_future_contacts(self):
        a,_=self.meshes();target=a.triangles+np.array([2.,0.,0.]);index=TriangleCrossIndex(target)
        self.assertEqual(list(index.batches(a.triangles)),[])
        pairs=np.concatenate(list(index.batches(a.triangles,padding_m=1.01)))
        np.testing.assert_array_equal(pairs,[[0,0]])
        with self.assertRaises(ValueError):list(index.batches(a.triangles,padding_m=-1))


if __name__ == '__main__':
    unittest.main()
