import unittest
import json
from types import SimpleNamespace
from unittest.mock import patch

import numpy as np

from geometry_tools.algorithms.section_course_socp import solve_weighted_section_socp


class BoundedWeightedCourseTests(unittest.TestCase):
    def test_nonfinite_failed_solver_evidence_is_saved_without_a_candidate(self):
        failed=SimpleNamespace(status='MaxTime',z=np.array([np.nan]),iterations=1,
            solve_time=.01,obj_val=np.nan,obj_val_dual=np.nan)
        with patch('geometry_tools.algorithms.section_course_socp.clarabel.DefaultSolver')as solver:
            solver.return_value.solve.return_value=failed
            d,proof=solve_weighted_section_socp(np.eye(1),[[1.]],[[1.,0.,0.]],[.001],maximum_motion_m=.003)
        self.assertIsNone(d);self.assertFalse(proof['feasible'])
        self.assertTrue(proof['nonfiniteSolverEvidence']);self.assertEqual(proof['status'],'MaxTime')
        self.assertIsNone(proof['solverPrimalObjective'])
        self.assertIsNone(proof['leadingHalfspaceDualRows'][0]['dual'])
        json.dumps(proof,allow_nan=False)

    def test_protected_infinite_radius_is_explicit_and_serializable(self):
        d,proof=solve_weighted_section_socp(np.diag([0.,1.]),[[1.,0.]],[[1.,0.,0.]],[.001],maximum_motion_m=.003)
        self.assertIsNone(d);self.assertTrue(proof['infiniteNecessaryRadius'])
        self.assertIsNone(proof['necessaryRadiusLowerBoundMeters'])
        json.dumps(proof,allow_nan=False)

    def test_elastic_route_reports_unsatisfied_original_rows_and_keeps_bone_hard(self):
        # A soft pair cannot both clear x=+/-1 mm. It must not turn that
        # diagnostic compromise into a claim that the original rows passed.
        d,proof=solve_weighted_section_socp(np.eye(1),[[1.],[1.],[1.]],
            [[1.,0.,0.],[-1.,0.,0.],[0.,1.,0.]],[.001,.001,.0005],
            maximum_motion_m=.003,elastic_row_groups=[0,0,-1])
        self.assertTrue(proof['feasible']);self.assertTrue(proof['proposalOnly'])
        self.assertFalse(proof['originalHalfspacesSatisfied'])
        self.assertGreater(proof['maximumUnrelaxedResidualMeters'],.00099)
        self.assertGreaterEqual(d[0,1],.0005-1e-10)
        self.assertLessEqual(np.linalg.norm(d[0]),.003+1e-10)

    def test_elastic_search_never_releases_an_incompatible_protected_bone_row(self):
        d,proof=solve_weighted_section_socp(np.diag([0.,1.]),[[1.,0.],[0.,1.]],
            [[1.,0.,0.],[1.,0.,0.]],[.001,.001],maximum_motion_m=.003,
            elastic_row_groups=[-1,0])
        self.assertIsNone(d);self.assertEqual(proof['status'],'ProtectedConstraintInfeasible')

    def test_exact_duplicate_halfspaces_keep_the_original_feasible_set(self):
        w=np.repeat([[1.,0.]],20,axis=0);n=np.repeat([[1.,0.,0.]],20,axis=0);b=np.full(20,.001)
        a,proof=solve_weighted_section_socp(np.eye(2),w,n,b,maximum_motion_m=.003)
        original,full=solve_weighted_section_socp(np.eye(2),w,n,b,maximum_motion_m=.003,deduplicate_halfspaces=False)
        self.assertTrue(proof['feasible']);self.assertTrue(full['feasible'])
        self.assertEqual(proof['exactDuplicateRowsRemoved'],19)
        np.testing.assert_allclose(a,original,atol=1e-9,rtol=0)
        self.assertLessEqual(np.max(b-np.einsum('ij,ij->i',n,w@a)),1e-10)
        self.assertTrue(all(x['row']==0 for x in proof['leadingHalfspaceDualRows']))

    def test_same_normal_with_different_bounds_is_not_deduplicated(self):
        a,proof=solve_weighted_section_socp(np.eye(1),[[1.],[1.]],[[1.,0.,0.],[1.,0.,0.]],
                                           [.001,.002],maximum_motion_m=.003)
        self.assertTrue(proof['feasible']);self.assertEqual(proof['exactDuplicateRowsRemoved'],0)
        self.assertAlmostEqual(a[0,0],.002,delta=1e-9)

    def test_explicit_sparse_precision_matches_the_same_dense_objective(self):
        H = np.array([[2., -1.], [-1., 2.]])
        K = np.linalg.inv(H)
        args = (K, [[.25, .75]], [[1., 0., 0.]], [.001])
        a, _ = solve_weighted_section_socp(*args, maximum_motion_m=.003)
        b, proof = solve_weighted_section_socp(*args, maximum_motion_m=.003, precision_matrix=H)
        self.assertTrue(proof['feasible'])
        np.testing.assert_allclose(a, b, atol=1e-10, rtol=0.)
        with self.assertRaisesRegex(ValueError, 'does not invert'):
            solve_weighted_section_socp(*args, maximum_motion_m=.003, precision_matrix=np.eye(2))

    def test_general_vector_guard_has_the_analytic_minimum_norm_solution(self):
        d, proof = solve_weighted_section_socp(np.eye(2), [[1., 0.]], [[1., 0., 0.]], [-.001],
            maximum_motion_m=.003, vector_rows=[[[0., 2., 0.], [0., 0., -1.]]], vector_offsets=[.001])
        self.assertTrue(proof['feasible'])
        np.testing.assert_allclose(d, [[0., .0004, 0.], [0., 0., -.0002]], atol=1e-8, rtol=0.)

    def test_shape_target_is_retained_when_inside_all_bounds(self):
        target = np.array([[.001, .002, 0.]])
        d, proof = solve_weighted_section_socp(np.eye(1), [[1.]], [[0., 0., 1.]], [-.001],
            maximum_motion_m=.003, target_displacement=target)
        self.assertTrue(proof['feasible'])
        np.testing.assert_allclose(d, target, atol=1e-8, rtol=0.)

    def test_ball_changes_the_coupled_metric_minimizer(self):
        # Without the radius the true minimizer is (0.27027, 3.24324) mm.
        # With a 3 mm ball it is exactly (1, 3) mm, satisfying the same
        # material-point constraint .25*d0 + .75*d1 >= 2.5 mm.
        K = np.diag([1., 4.])
        d, proof = solve_weighted_section_socp(K, [[.25, .75]], [[1., 0., 0.]],
            [.0025], maximum_motion_m=.003)
        self.assertTrue(proof['feasible'])
        np.testing.assert_allclose(d, [[.001, 0., 0.], [.003, 0., 0.]], atol=2e-8, rtol=0.)
        self.assertLessEqual(np.linalg.norm(d, axis=1).max(), .003 + 1e-10)

    def test_coupled_four_mm_requirement_is_infeasible_in_three_mm_balls(self):
        # A component-wise box would incorrectly permit this diagonal
        # motion. Euclidean balls cannot produce a weighted 4 mm offset.
        n = [[2. ** -.5, 2. ** -.5, 0.]]
        args = (np.eye(2), [[.5, .5]], n, [.004])
        d, proof = solve_weighted_section_socp(*args, maximum_motion_m=.003)
        self.assertIsNone(d)
        self.assertEqual(proof['status'], 'PrimalInfeasible')
        minimum, evidence = solve_weighted_section_socp(*args)
        self.assertTrue(evidence['feasible'])
        self.assertAlmostEqual(evidence['minimumRadiusMeters'], .004, delta=1e-9)
        self.assertAlmostEqual(evidence['minimumRadiusDualLowerBoundMeters'], .004, delta=1e-9)
        self.assertAlmostEqual(np.linalg.norm(minimum[0]), .004, delta=1e-9)

    def test_protected_station_is_literal_and_cannot_be_released(self):
        K = np.diag([0., 1.])
        d, proof = solve_weighted_section_socp(K, [[.5, .5]], [[0., 1., 0.]],
            [.001], maximum_motion_m=.003)
        self.assertTrue(proof['feasible'])
        np.testing.assert_array_equal(d[0], np.zeros(3))
        self.assertAlmostEqual(d[1, 1], .002, delta=1e-9)
        d, proof = solve_weighted_section_socp(K, [[1., 0.]], [[0., 1., 0.]],
            [.001], maximum_motion_m=.003)
        self.assertIsNone(d)
        self.assertEqual(proof['status'], 'ProtectedConstraintInfeasible')


if __name__ == '__main__':
    unittest.main()
