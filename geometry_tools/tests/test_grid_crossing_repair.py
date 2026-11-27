import unittest
import numpy as np
from geometry_tools.distance import triangle_pairs
from geometry_tools.algorithms.grid_crossing_repair import repair_crossing_rounding


def crossing_triangles():
    v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],
                [.2,.2,-1e-8],[.4,.2,1e-8],[.2,.4,1e-8]],np.float32).astype(float)
    f=np.arange(6).reshape(-1,3);p=v[f]
    n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
    return v,f,np.repeat(n[:,None],3,axis=1)


class GridCrossingRepairTests(unittest.TestCase):
    def test_crossing_is_removed_on_the_actual_float32_grid(self):
        v,f,w=crossing_triangles()
        self.assertEqual(triangle_pairs(v[f[:1]],v[f[1:]])['relation'][0],'crossing')
        r=repair_crossing_rounding(v,v,f,w,[[0,1]],origin=[0,0,0],locked=[0,1,2])
        self.assertTrue(r['evidence']['passed'])
        actual=r['positions']
        self.assertNotIn(triangle_pairs(actual[f[:1]],actual[f[1:]])['relation'][0],['crossing','coplanar-overlap'])
        np.testing.assert_array_equal(actual[:3],v[:3])
        np.testing.assert_array_equal(actual,actual.astype(np.float32).astype(float))
        self.assertLessEqual(np.linalg.norm(actual-v,axis=1).max(),5e-8)

    def test_locked_controls_do_not_move(self):
        v,f,w=crossing_triangles()
        r=repair_crossing_rounding(v,v,f,w,[[0,1]],origin=[0,0,0],locked=range(6))
        self.assertFalse(r['evidence']['passed'])
        np.testing.assert_array_equal(r['positions'],v)

    def test_source_bound_is_not_expanded_to_find_a_solution(self):
        v,f,w=crossing_triangles()
        r=repair_crossing_rounding(v,v,f,w,[[0,1]],origin=[0,0,0],maximum_displacement_m=1e-12)
        self.assertFalse(r['evidence']['passed'])
        np.testing.assert_array_equal(r['positions'],v)

    def test_combined_grid_search_keeps_the_same_source_bound(self):
        v,f,w=crossing_triangles()
        r=repair_crossing_rounding(v,v,f,w,[[0,1]],origin=[0,0,0],locked=[0,1,2],candidate_mode='bounded-combined-lattice')
        self.assertTrue(r['evidence']['passed'])
        self.assertLessEqual(np.linalg.norm(r['positions']-v,axis=1).max(),5e-8)


if __name__=='__main__':unittest.main()
