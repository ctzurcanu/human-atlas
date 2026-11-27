import unittest
import numpy as np
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.vertex_clearance_lp import repair_vertex_clearance


class VertexClearanceLPTests(unittest.TestCase):
    def fixture(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.2,.2,-1e-8],[.4,.2,1e-8],[.2,.4,1e-8]],np.float32).astype(float);f=np.arange(6).reshape(-1,3);p=v[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        return v,f,np.repeat(n[:,None],3,axis=1)

    def test_linear_clearance_preserves_fields_and_fixed_surface(self):
        v,f,w=self.fixture();guard=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        r=repair_vertex_clearance(v,v,f,w,origin=[0,0,0],global_guard=guard,maximum_displacement_m=5e-8,locked=[0,1,2])
        self.assertTrue(r['evidence']['passed']);self.assertEqual(len(guard.remaining),0)
        np.testing.assert_array_equal(r['positions'][:3],v[:3]);p=r['positions'][f];normal=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);self.assertTrue(np.all(np.einsum('fki,fi->fk',w,normal)>0))

    def test_locked_vertices_prevent_any_edit(self):
        v,f,w=self.fixture();guard=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        r=repair_vertex_clearance(v,v,f,w,origin=[0,0,0],global_guard=guard,maximum_displacement_m=5e-8,locked=range(6))
        self.assertFalse(r['evidence']['passed']);np.testing.assert_array_equal(r['positions'],v)

    def test_source_bound_is_not_expanded(self):
        v,f,w=self.fixture();guard=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        r=repair_vertex_clearance(v,v,f,w,origin=[0,0,0],global_guard=guard,maximum_displacement_m=1e-12,locked=[0,1,2])
        self.assertFalse(r['evidence']['passed']);np.testing.assert_array_equal(r['positions'],v)
