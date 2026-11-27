import unittest
import numpy as np
from geometry_tools.algorithms.grouped_float32_repair import repair_grouped_rounding
from geometry_tools.distance import triangle_pairs


class GroupedCrossingRepairTests(unittest.TestCase):
    def fixture(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.2,.2,-1e-8],[.4,.2,1e-8],[.2,.4,1e-8]],np.float32).astype(float)
        f=np.arange(6).reshape(-1,3);p=v[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        return v,f,np.repeat(n[:,None],3,axis=1)

    def test_explicit_rigid_triangle_translation_removes_a_crossing(self):
        v,f,w=self.fixture()
        r=repair_grouped_rounding(v,f,w,[[3,4,5]],origin=[0,0,0],stored_positions=v,maximum_displacement_m=5e-8,defect_face_pairs=[[0,1]],translations_m=[[0,0,3e-8]])
        self.assertTrue(r['evidence']['passed'])
        actual=r['positions'];np.testing.assert_array_equal(actual[:3],v[:3])
        self.assertLessEqual(np.ptp(actual[3:]-v[3:],axis=0).max(),1.5e-14)
        self.assertNotIn(triangle_pairs(actual[f[:1]],actual[f[1:]])['relation'][0],['crossing','coplanar-overlap'])
        self.assertEqual(r['evidence']['remaining_bad_faces'],0)

    def test_a_locked_member_protects_the_whole_group(self):
        v,f,w=self.fixture()
        r=repair_grouped_rounding(v,f,w,[[3,4,5]],origin=[0,0,0],stored_positions=v,maximum_displacement_m=5e-8,locked=[3],defect_face_pairs=[[0,1]],translations_m=[[0,0,3e-8]])
        self.assertFalse(r['evidence']['passed'])
        np.testing.assert_array_equal(r['positions'],v)

    def test_explicit_translation_cannot_expand_source_bound(self):
        v,f,w=self.fixture()
        r=repair_grouped_rounding(v,f,w,[[3,4,5]],origin=[0,0,0],stored_positions=v,maximum_displacement_m=5e-8,defect_face_pairs=[[0,1]],translations_m=[[0,0,1e-4]])
        self.assertFalse(r['evidence']['passed'])
        np.testing.assert_array_equal(r['positions'],v)

    def test_alias_cluster_bound_applies_to_every_original_source_member(self):
        v,f,w=self.fixture();v=np.r_[v,[v[3]+[0,0,-4e-8]]].astype(np.float32).astype(float)
        r=repair_grouped_rounding(v,f,w,[[3,4,5]],origin=[0,0,0],stored_positions=v,maximum_displacement_m=5e-8,
            defect_face_pairs=[[0,1]],translations_m=[[0,0,3e-8]],source_clusters={3:[3,6]},maximum_group_size=16)
        self.assertFalse(r['evidence']['passed']);np.testing.assert_array_equal(r['positions'],v)

    def test_larger_coherent_groups_require_an_explicit_bound(self):
        v,f,w=self.fixture()
        extra=np.array([p for k in range(4)for p in [[2+k,0,0],[2.25+k,0,0],[2+k,.25,0]]])
        v=np.r_[v,extra];f=np.r_[f,np.arange(6,18).reshape(-1,3)];w=np.r_[w,np.tile([[[0.,0.,1.]]*3],(4,1,1))]
        args=dict(origin=[0,0,0],stored_positions=v,maximum_displacement_m=5e-8,defect_face_pairs=[[0,1]],translations_m=[[0,0,3e-8]])
        group=[3,4,5,*range(6,18)]
        with self.assertRaises(ValueError):repair_grouped_rounding(v,f,w,[group],**args)
        r=repair_grouped_rounding(v,f,w,[group],maximum_group_size=16,**args)
        self.assertTrue(r['evidence']['passed'])

    def test_unused_control_cannot_enter_a_moving_group(self):
        v,f,w=self.fixture();v=np.r_[v,[[0.,0.,2.]]]
        with self.assertRaises(ValueError):repair_grouped_rounding(v,f,w,[[3,4,5,6]],origin=[0,0,0],maximum_displacement_m=5e-8,translations_m=[[0,0,3e-8]])


if __name__=='__main__':unittest.main()
