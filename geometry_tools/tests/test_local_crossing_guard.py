import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH,box_distance_squared
from geometry_tools.distance import triangle_pairs
from geometry_tools.algorithms.local_crossing_guard import LocalCrossingGuard
from geometry_tools.algorithms.grid_crossing_repair import repair_crossing_rounding
from geometry_tools.algorithms.grouped_float32_repair import repair_grouped_rounding
from geometry_tools.algorithms.local_crossing_guard import CrossingBudgetExceeded


def scene():
    v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],
        [.2,.2,-1e-8],[.4,.2,1e-8],[.2,.4,1e-8],
        [0.,0.,-2e-8],[1.,0.,-2e-8],[0.,1.,-2e-8]])
    return v,np.arange(9).reshape(-1,3)


class LocalCrossingGuardTests(unittest.TestCase):
    def test_box_queries_match_exhaustive_boxes_with_and_without_margin(self):
        rng=np.random.default_rng(918);v=rng.random((150,3));f=np.arange(150).reshape(-1,3);mesh=Mesh(v,f,'test');tree=BVH(mesh,leaf_size=4)
        for margin in [0.,.03]:
            for _ in range(8):
                lo=rng.random(3);hi=lo+.04;t=mesh.triangles
                expected=np.flatnonzero(box_distance_squared(lo,hi,t.min(1),t.max(1))<=margin*margin)
                np.testing.assert_array_equal(tree.query_bounds(lo,hi,distance_m=margin),expected)

    def test_fixing_a_known_crossing_cannot_create_one_with_another_triangle(self):
        v,f=scene();g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        ids=np.array([3,4,5]);bad=v[ids]+[0,0,-2.5e-8]
        self.assertIsNone(g.check_move(ids,bad),'repairing the old pair would cross the second plane')
        self.assertEqual(len(g.remaining),1);np.testing.assert_array_equal(g.positions,v)
        good=g.check_move(ids,v[ids]+[0,0,3e-8]);self.assertIsNotNone(good);g.commit(good)
        self.assertEqual(len(g.remaining),0)
        pairs=np.array([[0,1],[0,2],[1,2]])
        self.assertFalse(np.isin(triangle_pairs(g.positions[f[pairs[:,0]]],g.positions[f[pairs[:,1]]])['relation'],['crossing','coplanar-overlap']).any())
        with self.assertRaises(ValueError):g.commit(good)

    def test_source_bounds_and_literal_index_are_enforced(self):
        v,f=scene();g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        with self.assertRaises(ValueError):g.check_move([3],v[[3]]+[0,0,2e-7])
        with self.assertRaises(ValueError):LocalCrossingGuard(v,f,[[0,2]],maximum_displacement_m=1e-7)
        other=Mesh(v+[.01,0,0],f,'test')
        with self.assertRaises(ValueError):LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7,bvh=BVH(other))

    def test_grid_repair_runs_with_global_protection(self):
        v,f=scene();v=v.astype(np.float32).astype(float);p=v[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        r=repair_crossing_rounding(v,v,f,np.repeat(n[:,None],3,axis=1),[[0,1]],origin=[0,0,0],locked=[0,1,2,6,7,8],global_guard=g)
        self.assertTrue(r['evidence']['passed']);self.assertEqual(len(g.remaining),0)
        np.testing.assert_array_equal(r['positions'],g.positions)
        self.assertFalse(r['evidence']['global_crossing_guard']['new_global_defects_permitted'])

    def test_deforming_grid_repair_emits_bounded_transported_outward_words(self):
        v,f=scene();v=v.astype(np.float32).astype(float);p=v[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        words=np.rint(np.repeat(n[:,None],3,axis=1)*32767).astype(np.int16)
        g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        r=repair_crossing_rounding(v,v,f,words,[[0,1]],origin=[0,0,0],locked=[0,1,2,6,7,8],global_guard=g,maximum_facet_rotation_degrees=45)
        self.assertTrue(r['evidence']['passed']);self.assertIsNotNone(r['evidence']['normal_transport'])
        p=r['positions'][f];normal=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
        self.assertTrue(np.all(np.einsum('fki,fi->fk',r['corner_words'].astype(float),normal)>0))
        np.testing.assert_array_equal(r['positions'][[0,1,2,6,7,8]],v[[0,1,2,6,7,8]])

    def test_group_search_skips_a_closer_translation_that_creates_a_new_crossing(self):
        v,f=scene();v=v.astype(np.float32).astype(float);p=v[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        r=repair_grouped_rounding(v,f,np.repeat(n[:,None],3,axis=1),[[3,4,5]],origin=[0,0,0],stored_positions=v,
            maximum_displacement_m=5e-8,defect_face_pairs=[[0,1]],translations_m=[[0,0,-2.5e-8],[0,0,3e-8]],global_guard=g)
        self.assertTrue(r['evidence']['passed']);self.assertEqual(len(g.remaining),0)
        self.assertGreater(r['positions'][3,2],0)
        self.assertEqual(g.rejected_new,1)

    def test_pair_budget_rejects_without_changing_the_geometry(self):
        v,f=scene();g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7,maximum_pair_tests=1)
        with self.assertRaises(CrossingBudgetExceeded):g.check_move([3,4,5],v[[3,4,5]]+[0,0,3e-8])
        np.testing.assert_array_equal(g.positions,v)

    def test_retriangulation_checks_all_neighbors_and_keeps_stable_face_ids(self):
        v,f=scene();v=np.r_[v,v[[3,4,5]]+[0,0,-2.5e-8],v[[3,4,5]]+[0,0,3e-8]]
        g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        self.assertIsNone(g.check_faces([1],[[9,10,11]],[True]))
        token=g.check_faces([1],[[12,13,14]],[True]);self.assertIsNotNone(token);g.commit(token)
        self.assertEqual(len(g.remaining),0);np.testing.assert_array_equal(g.faces[1],[12,13,14])
        np.testing.assert_array_equal(g.positions,v)
        with self.assertRaises(ValueError):g.check_move([3],v[[3]])

    def test_deleted_faces_are_excluded_and_cannot_be_reactivated(self):
        v,f=scene();g=LocalCrossingGuard(v,f,[[0,1]],maximum_displacement_m=1e-7)
        token=g.check_faces([1],f[[1]],[False]);g.commit(token)
        self.assertEqual(len(g.remaining),0)
        with self.assertRaises(ValueError):g.check_faces([1],f[[1]],[True])

    def test_disjoint_actual_boxes_are_excluded_after_conservative_source_query(self):
        first=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        v=np.r_[first,first+[0,0,1e-17]];f=np.arange(6).reshape(-1,3)
        g=LocalCrossingGuard(v,f,[],maximum_displacement_m=1e-7)
        proposal=g.check_move([3,4,5],v[[3,4,5]]+[1e-9,0,0]);self.assertIsNotNone(proposal)
        g.commit(proposal);self.assertEqual(len(g.remaining),0);self.assertEqual(g.tests,0)


if __name__=='__main__':unittest.main()
