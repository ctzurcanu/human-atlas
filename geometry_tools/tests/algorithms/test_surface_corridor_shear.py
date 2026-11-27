import unittest
import numpy as np
from geometry_tools.algorithms.surface_corridor_shear import fit_vertex_corridor,fit_z_corridor,basis


class LocalCorridorTests(unittest.TestCase):
    @staticmethod
    def finite_box():
        v=np.array([[0.,0.,0.],[1.,0.,0.],[1.,1.,0.],[0.,1.,0.],
                    [0.,0.,1.],[1.,0.,1.],[1.,1.,1.],[0.,1.,1.]])
        f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],
                    [0,1,5],[0,5,4],[1,2,6],[1,6,5],
                    [2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        return v,f

    def test_declared_adjacent_ceiling_allows_separation_but_blocks_penetration(self):
        v,f=self.finite_box();ceiling=v[[[4,6,5],[4,7,6]]].copy()
        moving_ceiling=ceiling.copy();moving_ceiling[:,:,2]=.8
        free=fit_vertex_corridor(v,f,[],[('new-ceiling',moving_ceiling)],
            preserve_peers=[('contact',ceiling)],declared_contact_peers=['contact'])
        self.assertEqual(free['status'],'Solved')
        np.testing.assert_array_equal(free['positions'][:,:2],v[:,:2])
        np.testing.assert_allclose(free['positions'][4:,2],.8,atol=1e-8,rtol=0)
        np.testing.assert_allclose(free['positions'][4:,2]-free['positions'][:4,2],1.,atol=1e-9,rtol=0)
        self.assertGreater(free['counts']['contact-declared-zero-gap-contacts'],0)
        floor=v[[[0,1,2],[0,2,3]]].copy();floor[:,:,2]=.2
        blocked=fit_vertex_corridor(v,f,[('raised-floor',floor)],[],
            preserve_peers=[('contact',ceiling)],declared_contact_peers=['contact'])
        self.assertNotEqual(blocked['status'],'Solved')
        self.assertFalse(blocked['independentLinearFeasibility'])

    def test_declared_adjacent_floor_retains_the_opposite_material_side(self):
        v,f=self.finite_box();floor=v[[[0,1,2],[0,2,3]]].copy()
        ceiling=v[[[4,6,5],[4,7,6]]].copy();ceiling[:,:,2]=.8
        blocked=fit_vertex_corridor(v,f,[],[('lowered-ceiling',ceiling)],
            preserve_peers=[('contact',floor)],declared_contact_peers=['contact'])
        self.assertNotEqual(blocked['status'],'Solved')
        self.assertFalse(blocked['independentLinearFeasibility'])

    def test_contact_must_be_declared_and_oppositely_oriented(self):
        v,f=self.finite_box();contact=v[[[4,6,5],[4,7,6]]].copy()
        default=fit_vertex_corridor(v,f,[],[],preserve_peers=[('contact',contact)])
        self.assertEqual(default['status'],'OriginalPeerOrderAmbiguous')
        wrong=fit_vertex_corridor(v,f,[],[],preserve_peers=[('contact',contact[:,[0,2,1]])],
            declared_contact_peers=['contact'])
        self.assertEqual(wrong['status'],'OriginalPeerOrderAmbiguous')
        crossing=contact.copy();crossing[0,0,2]=.5;crossing[0,1,2]=1.5
        wrong=fit_vertex_corridor(v,f,[],[],preserve_peers=[('contact',crossing)],
            declared_contact_peers=['contact'])
        self.assertEqual(wrong['status'],'OriginalPeerOrderAmbiguous')
        with self.assertRaises(ValueError):
            fit_vertex_corridor(v,f,[],[],declared_contact_peers=['absent'])

    def test_original_peer_side_cannot_be_crossed_to_satisfy_another_contact(self):
        v = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.]])
        f = np.array([[0, 1, 2], [0, 2, 1]])
        lower = v[[[0, 1, 2]]].copy(); lower[:, :, 2] = 1.
        upper = v[[[0, 2, 1]]].copy(); upper[:, :, 2] = 2.
        obstacle = v[[[0, 1, 2]]].copy(); obstacle[:, :, 2] = .5
        free = fit_vertex_corridor(v, f, [('lower', lower)], [('upper', upper)], guard=.01)
        self.assertEqual(free['status'], 'Solved')
        protected = fit_vertex_corridor(v, f, [('lower', lower)], [('upper', upper)], guard=.01,
                                       preserve_peers=[('obstacle', obstacle)])
        self.assertNotEqual(protected['status'], 'Solved')
        self.assertFalse(protected['independentLinearFeasibility'])

    def test_local_scope_excludes_only_literal_pinned_triangles(self):
        # The first reference patch has an inherited violation. Its controls
        # are fixed; the second patch must actually move below its upper peer.
        v = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.],
                      [2., 0., 2.], [3., 0., 2.], [2., 1., 2.]])
        f = np.array([[0, 1, 2], [3, 4, 5]])
        pinned = np.array([True, True, True, False, False, False])
        lower = v[[[0, 1, 2]]].copy(); lower[:, :, 2] = 1.
        upper = v[[[3, 5, 4]]].copy(); upper[:, :, 2] = 1.
        doubled = np.r_[f, f[:, [0, 2, 1]]]
        strict = fit_vertex_corridor(v, doubled, [('lower', lower)], [('upper', upper)],
                                     pinned=pinned, guard=.01)
        self.assertEqual(strict['status'], 'InfeasibleFixedTriangle')
        local = fit_vertex_corridor(v, doubled, [('lower', lower)], [('upper', upper)],
                                    pinned=pinned, guard=.01, moving_faces_only=True)
        self.assertEqual(local['status'], 'Solved')
        np.testing.assert_array_equal(local['positions'][:3], v[:3])
        np.testing.assert_allclose(local['positions'][3:, :2], v[3:, :2])
        self.assertTrue(np.all(local['positions'][3:, 2] <= .99+1e-9))
        self.assertGreater(local['counts']['upper'], 0)

    def test_local_spline_support_preserves_distant_controls_and_direction(self):
        v=np.array([[0.,0.,1.],[1.,0.,1.],[0.,1.,1.],[10.,10.,4.],[11.,10.,4.],[10.,11.,4.]])
        f=np.array([[0,1,2],[3,4,5]]);ceiling=v[[[0,2,1]]].copy();ceiling[:,:,2]=.8
        result=fit_z_corridor(v,f,[],[('ceiling',ceiling)],spacing=.5,pinned=np.array([False]*3+[True]*3),
            posterior_only=False,knot_bounds=[[-1.,-1.],[2.,2.]],nonpositive_vertices=True,coefficient_bound=.5)
        self.assertEqual(result['status'],'Solved')
        np.testing.assert_array_equal(basis(v[3:],result['knots'],.5),0.)
        np.testing.assert_array_equal(result['positions'][3:],v[3:])
        self.assertTrue(np.all(result['displacements']<=1e-8))
        self.assertTrue(np.all(result['positions'][:3,2]<=.8+1e-8))
        with self.assertRaises(ValueError):
            fit_z_corridor(v,f,[],[('ceiling',ceiling)],spacing=.5,knot_bounds=[[1.,1.],[0.,0.]])

    def test_weighted_spline_compresses_only_the_declared_surface(self):
        v,f=self.finite_box();ceiling=v[[[4,6,5],[4,7,6]]].copy();ceiling[:,:,2]=.8
        weights=np.r_[np.zeros(4),np.ones(4)]
        result=fit_z_corridor(v,f,[],[('ceiling',ceiling)],spacing=.5,
            posterior_only=False,nonpositive_vertices=True,coefficient_bound=.5,displacement_weights=weights,solver_variable_scale=.1)
        self.assertEqual(result['status'],'Solved')
        np.testing.assert_array_equal(result['positions'][:4],v[:4])
        self.assertIsNone(result['continuousJacobianDeterminant'])
        self.assertIsNone(result['gradient'])
        self.assertTrue(result['weightedFieldDerivativeRequired'])
        self.assertEqual(result['solverVariableScale'],.1)
        np.testing.assert_array_equal(result['positions'][:,:2],v[:,:2])
        self.assertTrue(np.all(result['positions'][4:,2]<=.8+1e-8))
        self.assertTrue(np.all(result['positions'][4:,2]>.5))
        for invalid in [weights[:-1],np.r_[weights[:-1],np.nan],weights+1.]:
            with self.assertRaises(ValueError):
                fit_z_corridor(v,f,[],[('ceiling',ceiling)],spacing=.5,displacement_weights=invalid)

    def test_local_spline_refit_omits_only_facets_with_zero_response(self):
        v=np.array([[0.,0.,1.],[1.,0.,1.],[0.,1.,1.],[10.,10.,4.],[11.,10.,4.],[10.,11.,4.]])
        f=np.array([[0,1,2],[3,4,5]])
        ceiling=v[[[0,2,1]]].copy();ceiling[:,:,2]=.8
        fixed=v[[[3,5,4]]].copy();fixed[:,:,2]=3.
        arguments=dict(spacing=.5,posterior_only=False,knot_bounds=[[-1.,-1.],[2.,2.]],coefficient_bound=.5)
        strict=fit_z_corridor(v,f,[],[('ceiling',ceiling),('inherited-fixed-reference',fixed)],**arguments)
        self.assertEqual(strict['status'],'InfeasibleFixedRow')
        local=fit_z_corridor(v,f,[],[('ceiling',ceiling),('inherited-fixed-reference',fixed)],moving_faces_only=True,**arguments)
        self.assertEqual(local['status'],'Solved')
        np.testing.assert_array_equal(local['positions'][3:],v[3:])
        self.assertTrue(np.all(local['positions'][:3,2]<=.8+1e-8))

    def test_peer_clearance_keeps_a_literal_floor_while_moving_the_ceiling(self):
        v,f=self.finite_box();ceiling=v[[[4,6,5],[4,7,6]]].copy();ceiling[:,:,2]=.8
        floor=v[[[0,1,2],[0,2,3]]].copy();weights=np.r_[np.zeros(4),np.ones(4)]
        kw=dict(spacing=.5,posterior_only=False,coefficient_bound=.5,displacement_weights=weights,guard=.01)
        strict=fit_z_corridor(v,f,[('literal-floor',floor)],[('ceiling',ceiling)],**kw)
        self.assertEqual(strict['status'],'InfeasibleFixedRow')
        self.assertEqual(strict['fixedConstraintExamples'][0]['peer'],'literal-floor')
        fit=fit_z_corridor(v,f,[('literal-floor',floor)],[('ceiling',ceiling)],peer_guards={'literal-floor':0.},**kw)
        self.assertEqual(fit['status'],'Solved');np.testing.assert_array_equal(fit['positions'][:4],v[:4])
        self.assertTrue(np.all(fit['positions'][4:,2]<=.79+1e-8))
        with self.assertRaises(ValueError):fit_z_corridor(v,f,[],[('ceiling',ceiling)],spacing=.5,peer_guards={'absent':0.})

    def test_compound_peer_uses_declared_real_pairs_instead_of_one_global_side(self):
        v=np.array([[0.,0.,1.],[1.,0.,1.],[0.,1.,1.],[0.,0.,3.],[1.,0.,3.],[0.,1.,3.]])
        f=np.array([[0,1,2],[3,4,5]]);peer=v[[[0,2,1],[3,5,4]]].copy();peer[0,:,2]=0.;peer[1,:,2]=2.
        kw=dict(spacing=.5,posterior_only=False,coefficient_bound=1.5,displacement_weights=np.r_[np.zeros(3),np.ones(3)])
        strict=fit_z_corridor(v,f,[],[('folded-peer',peer)],**kw)
        self.assertEqual(strict['status'],'InfeasibleFixedRow')
        local=fit_z_corridor(v,f,[],[('folded-peer',peer)],triangle_pairs={'folded-peer':np.array([[1,1]])},**kw)
        self.assertEqual(local['status'],'Solved');np.testing.assert_array_equal(local['positions'][:3],v[:3])
        self.assertTrue(np.all(local['positions'][3:,2]<=2.+1e-8))
        with self.assertRaises(ValueError):fit_z_corridor(v,f,[],[('folded-peer',peer)],spacing=.5,triangle_pairs={'folded-peer':np.array([[99,1]])})


if __name__ == '__main__':
    unittest.main()
