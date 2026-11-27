import unittest
import numpy as np
from scipy import sparse
from geometry_tools.algorithms.sequential_original_shape_sections import solve_sequential_original_shape


class SequentialSourceShape(unittest.TestCase):
    def test_additional_exact_contact_constraint_cannot_be_skipped(self):
        def shape(d):
            return sparse.csr_matrix((0,3)),np.array([]),dict(allOriginalFacetsChecked=1,violatingOriginalFacets=0)
        def contacts(d):
            bad=d[0,0]>.0001+1e-10
            return (sparse.csr_matrix([[-1.,0.,0.]]) if bad else sparse.csr_matrix((0,3)),np.array([-.0001]) if bad else np.array([]),dict(complete=True,violatingExactConstraints=int(bad)))
        d,proof=solve_sequential_original_shape(sparse.eye(1),sparse.csr_matrix((0,3)),np.array([]),sparse.eye(3),np.array([.0008,0.,0.]),exact_source_shape_rows=shape,additional_exact_rows=contacts,maximum_motion_m=.001,maximum_stages=4,wall_seconds=8.)
        self.assertIsNotNone(d);self.assertLessEqual(d[0,0],.0001+1e-10)
        self.assertGreater(len(proof['stages']),1)
    def test_source_limit_and_preexisting_context_both_replayed(self):
        def source(d):
            bad = d[0, 0] > .0001+1e-10
            return (sparse.csr_matrix([[-1., 0., 0.]]) if bad else sparse.csr_matrix((0, 3)),
                    np.array([-.0001]) if bad else np.array([]),
                    dict(allOriginalFacetsChecked=1, violatingOriginalFacets=int(bad)))
        A = sparse.csr_matrix([[0., -1., 0.]])
        d, proof = solve_sequential_original_shape(sparse.eye(1), A, np.array([-.0002]),
            sparse.eye(3), np.array([.0008, .0007, 0.]), exact_source_shape_rows=source,
            maximum_motion_m=.001, maximum_stages=4, wall_seconds=8.)
        self.assertIsNotNone(d)
        self.assertLessEqual(d[0, 0], .0001+1e-10)
        self.assertLessEqual(d[0, 1], .0002+1e-10)
        self.assertTrue(proof['fullOriginalSourceShapeChecked'])

    def test_exhaustion_never_emits_unchecked_shape(self):
        def source(d):
            return sparse.csr_matrix([[-1., 0., 0.]]), np.array([-.0001]), dict(
                allOriginalFacetsChecked=1, violatingOriginalFacets=1)
        d, proof = solve_sequential_original_shape(sparse.eye(1), sparse.csr_matrix((0, 3)), np.array([]),
            sparse.eye(3), np.array([.0008, 0., 0.]), exact_source_shape_rows=source,
            maximum_motion_m=.001, maximum_stages=1, wall_seconds=8.)
        self.assertIsNone(d)
        self.assertFalse(proof['feasible'])

    def test_completed_constraints_resume_without_accepting_incomplete_point(self):
        saved = {}
        def source(d):
            bad = d[0, 0] > .0001+1e-10
            return (sparse.csr_matrix([[-1., 0., 0.]]) if bad else sparse.csr_matrix((0, 3)),
                    np.array([-.0001]) if bad else np.array([]),
                    dict(allOriginalFacetsChecked=1, violatingOriginalFacets=int(bad)))
        def checkpoint(stage, d, matrix, rhs, seed, proof):
            saved.update(matrix=matrix.copy(), rhs=rhs.copy(), seed=seed.copy())
        args = (sparse.eye(1), sparse.csr_matrix((0, 3)), np.array([]), sparse.eye(3), np.array([.0008, 0., 0.]))
        d, proof = solve_sequential_original_shape(*args, exact_source_shape_rows=source,
            maximum_motion_m=.001, maximum_stages=1, wall_seconds=8., state_checkpoint=checkpoint)
        self.assertIsNone(d)
        self.assertFalse(proof['feasible'])
        resumed, proof = solve_sequential_original_shape(*args, exact_source_shape_rows=source,
            maximum_motion_m=.001, maximum_stages=2, wall_seconds=8.,
            initial_local_rows=(saved['matrix'], saved['rhs']), initial_active_ids=saved['seed'])
        self.assertIsNotNone(resumed)
        self.assertLessEqual(resumed[0, 0], .0001+1e-10)
        self.assertTrue(proof['fullOriginalSourceShapeChecked'])


if __name__ == '__main__':
    unittest.main()
