import unittest
import numpy as np
from geometry_tools.algorithms.grouped_float32_repair import repair_grouped_rounding


class GroupedFloat32RepairTests(unittest.TestCase):
    def fixture(self):
        v=np.array([[1.,0.,0.],[1.+4e-8,0.,0.],[1.,1.,0.],[1.,0.,1.]])
        f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
        p=v[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
        return v,f,np.repeat(np.rint(n*32767).astype(np.int16)[:,None],3,axis=1)

    def test_rigid_pair_opens_a_collapsed_tetrahedron_without_changing_fields(self):
        v,f,w=self.fixture();old_f=f.copy();old_w=w.copy();initial=v.astype(np.float32).astype(float)
        result=repair_grouped_rounding(v,f,w,[(0,2)],origin=[0,0,0],maximum_displacement_m=3e-7,wall_seconds=2)
        e=result['evidence'];self.assertTrue(e['passed']);self.assertLessEqual(e['maximum_displacement_m'],3e-7)
        shift=result['positions'][[0,2]]-initial[[0,2]]
        np.testing.assert_allclose(shift[0],shift[1],atol=1e-14)
        np.testing.assert_array_equal(f,old_f);np.testing.assert_array_equal(w,old_w)
        self.assertTrue(all(a>b for a,b in zip([e['initial_bad_faces']]+[m['remaining_bad_faces']for m in e['moves']], [m['remaining_bad_faces']for m in e['moves']])))

    def test_locked_group_is_not_moved(self):
        v,f,w=self.fixture();r=repair_grouped_rounding(v,f,w,[(0,2)],origin=[0,0,0],maximum_displacement_m=3e-7,locked=[2],wall_seconds=2)
        self.assertFalse(r['evidence']['passed']);self.assertEqual(r['evidence']['moves'],[])
        np.testing.assert_array_equal(r['positions'],v.astype(np.float32).astype(float))

    def test_zero_evaluation_budget_is_not_extended(self):
        v,f,w=self.fixture();r=repair_grouped_rounding(v,f,w,[(0,2)],origin=[0,0,0],maximum_displacement_m=3e-7,max_evaluations=0,wall_seconds=2)
        self.assertFalse(r['evidence']['passed']);self.assertEqual(r['evidence']['evaluations'],0)
        self.assertEqual(r['evidence']['budget_or_stop_reason'],'evaluation budget')

    def test_out_of_bounds_initial_candidate_is_rejected(self):
        v,f,w=self.fixture();q=v.astype(np.float32).astype(float);q[0,0]+=1.
        with self.assertRaisesRegex(ValueError,'exceeds'):
            repair_grouped_rounding(v,f,w,[(0,2)],origin=[0,0,0],maximum_displacement_m=3e-7,stored_positions=q)
