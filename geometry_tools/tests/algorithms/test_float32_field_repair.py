import unittest
import numpy as np
from geometry_tools.algorithms.float32_field_repair import repair_outward_rounding


class Float32FieldRepairTests(unittest.TestCase):
    def fixture(self):
        v=np.array([[1.,0.,0.],[1.+4e-8,0.,0.],[1.,1.,0.],[1.,0.,1.]])
        f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
        p=v[f];n=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])
        n/=np.linalg.norm(n,axis=1)[:,None]
        words=np.repeat(np.rint(n*32767).astype(np.int16)[:,None],3,axis=1)
        return v,f,words

    def test_collapsed_rounding_is_repaired_without_changing_faces_or_fields(self):
        v,f,w=self.fixture();old_f=f.copy();old_w=w.copy()
        result=repair_outward_rounding(v,f,w,maximum_displacement_m=3e-7,wall_seconds=2)
        self.assertTrue(result['evidence']['passed'])
        self.assertGreater(result['evidence']['initial_bad_faces'],0)
        self.assertEqual(result['evidence']['remaining_bad_faces'],0)
        self.assertTrue(np.array_equal(f,old_f));self.assertTrue(np.array_equal(w,old_w))
        self.assertLessEqual(result['evidence']['maximum_displacement_m'],3e-7)
        self.assertTrue(np.array_equal(result['relative_positions'].astype(float),result['positions']))

    def test_locked_source_is_not_moved_to_hide_a_failure(self):
        v,f,w=self.fixture()
        result=repair_outward_rounding(v,f,w,maximum_displacement_m=3e-7,locked=range(4),wall_seconds=2)
        self.assertFalse(result['evidence']['passed'])
        self.assertEqual(result['evidence']['moves'],[])
        self.assertTrue(np.array_equal(result['positions'],v.astype(np.float32).astype(float)))

    def test_explicit_evaluation_budget_is_reported(self):
        v,f,w=self.fixture()
        result=repair_outward_rounding(v,f,w,maximum_displacement_m=3e-7,max_evaluations=0,wall_seconds=2)
        self.assertFalse(result['evidence']['passed'])
        self.assertEqual(result['evidence']['budget_or_stop_reason'],'evaluation budget')
        self.assertEqual(result['evidence']['evaluations'],0)

    def test_first_improvement_keeps_the_same_safety_checks(self):
        v,f,w=self.fixture();before_f=f.copy();before_w=w.copy()
        result=repair_outward_rounding(v,f,w,maximum_displacement_m=3e-7,
            selection='first-improving',wall_seconds=2)
        evidence=result['evidence']
        self.assertTrue(evidence['passed'])
        self.assertEqual(evidence['search_selection'],'first-improving')
        self.assertLessEqual(evidence['maximum_displacement_m'],3e-7)
        counts=[evidence['initial_bad_faces']]+[m['remaining_bad_faces']for m in evidence['moves']]
        self.assertTrue(all(a>=b for a,b in zip(counts,counts[1:])))
        np.testing.assert_array_equal(f,before_f);np.testing.assert_array_equal(w,before_w)

    def test_unknown_search_selection_is_rejected(self):
        v,f,w=self.fixture()
        with self.assertRaisesRegex(ValueError,'selection'):
            repair_outward_rounding(v,f,w,maximum_displacement_m=3e-7,selection='unbounded')

    def test_bounded_lattice_moves_strictly_reduce_bad_face_count(self):
        v,f,w=self.fixture()
        result=repair_outward_rounding(v,f,w,maximum_displacement_m=3e-7,
            selection='first-improving',candidate_mode='bounded-global-lattice',
            require_bad_face_reduction=True,wall_seconds=2)
        evidence=result['evidence']
        self.assertTrue(evidence['passed'])
        counts=[evidence['initial_bad_faces']]+[m['remaining_bad_faces']for m in evidence['moves']]
        self.assertTrue(all(a>b for a,b in zip(counts,counts[1:])))
        self.assertLessEqual(evidence['maximum_displacement_m'],3e-7)
        np.testing.assert_array_equal(result['positions'],result['relative_positions'].astype(float))

    def test_lattice_radius_cannot_expand_without_a_finite_limit(self):
        v,f,w=self.fixture()
        with self.assertRaisesRegex(ValueError,'seven global quanta'):
            repair_outward_rounding(v,f,w,maximum_displacement_m=1.,
                candidate_mode='bounded-global-lattice')


if __name__=='__main__':unittest.main()
