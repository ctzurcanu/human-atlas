import unittest
import numpy as np
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.tree_section_corridor import tree_section_corridor,source_span_hits
from geometry_tools.tests.test_sparse_group_clearance import cube


class SourceCorridorTests(unittest.TestCase):
    def test_protected_mouth_is_not_released_to_make_a_grid_route_pass(self):
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]));skin=SurfaceIndex(cube([-.02,-.02,-.02],[.03,.03,.03]))
        v=np.array([[.0099,.004,.004],[.0099,.006,.004],[.0099,.005,.006]])
        f=np.array([[0,1,2]]);labels=np.zeros(3,int)
        _,free=tree_section_corridor(v,f,labels,bone,skin)
        q,fixed=tree_section_corridor(v,f,labels,bone,skin,protected_raw=[0,1,2])
        self.assertTrue(free['complete'],free);self.assertFalse(fixed['complete'],fixed)
        np.testing.assert_array_equal(q,v)

    def test_fine_local_grid_uses_saved_source_closure_without_dense_global_grid(self):
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]));skin=SurfaceIndex(cube([-.02,-.02,-.02],[.03,.03,.03]))
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006]])
        q,proof=tree_section_corridor(v,np.array([[0,1,2]]),np.zeros(3,int),bone,skin,
                    guide_vertices=v+[.001,0.,0.],grid_spacing_m=.00025,local_radius_m=.00075)
        self.assertTrue(proof['complete'],proof)
        self.assertTrue(proof['savedRouteGuideUsed']);self.assertEqual(proof['gridStates'],123)
        np.testing.assert_allclose(q-v,np.broadcast_to(q[0]-v[0],q.shape),atol=1e-15)

    def test_span_screen_finds_crossing_between_clear_vertices(self):
        surface=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        triangles=np.array([[[-.01,.005,.005],[.02,.005,.005],[-.01,.006,.006]],
                            [[.012,.004,.004],[.012,.006,.004],[.012,.005,.006]]])
        np.testing.assert_array_equal(source_span_hits(triangles,surface),[True,False])

    def test_closed_recipient_context_rejects_contained_states_and_preserves_source(self):
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        skin=SurfaceIndex(cube([-.02,-.02,-.02],[.03,.03,.03]))
        recipient=SurfaceIndex(cube([-.02,-.02,-.02],[.04,.04,.04]))
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006]])
        q,proof=tree_section_corridor(v,np.array([[0,1,2]]),np.zeros(3,int),bone,skin,
                    recipient_surfaces=[(recipient,True,'fixed tendon diagnostic')])
        self.assertFalse(proof['complete'],proof)
        self.assertEqual(proof['status'],'NoSampledRigidGroupCorridor')
        np.testing.assert_array_equal(q,v)

    def test_source_cycle_moves_as_one_rigid_closure(self):
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        skin=SurfaceIndex(cube([-.02,-.02,-.02],[.03,.03,.03]))
        v=np.array([[.0101,.004,.004],[.0101,.006,.004],[.0101,.005,.006]])
        q,proof=tree_section_corridor(v,np.array([[0,1,2]]),np.array([0,1,2]),bone,skin)
        self.assertTrue(proof['complete'],proof)
        self.assertEqual(proof['sourceGroups'],1)
        self.assertEqual(proof['originalSourceGroups'],3)
        np.testing.assert_allclose(q-v,np.broadcast_to(q[0]-v[0],q.shape),atol=1e-15)

    def test_no_sampled_corridor_preserves_original(self):
        surface=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        v=np.array([[.0101,.004,.004],[.0101,.006,.004],[.0101,.005,.006]])
        q,proof=tree_section_corridor(v,np.array([[0,1,2]]),np.array([0,0,0]),surface,surface)
        self.assertFalse(proof['complete'])
        np.testing.assert_array_equal(q,v)


if __name__=='__main__':unittest.main()
