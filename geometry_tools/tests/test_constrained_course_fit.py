import unittest
import numpy as np
from geometry_tools.algorithms.constrained_course_fit import fit_course


class ConstrainedCourseFitTest(unittest.TestCase):
    def fixture(self,offset):
        c=np.array([[0.,0.,0.],[0.,.01,0.],[0.,.02,0.]])
        h=[dict(coefficients={1:[1.,0.]},minimum_m=offset),dict(coefficients={1:[-1.,0.]},minimum_m=-offset)]
        return c,h

    def test_current_arc_iteration_resolves_a_false_source_parameter_rejection(self):
        c,h=self.fixture(.00464);original=c.copy();saved=[]
        delta,r=fit_course(c,[[0,1],[1,2]],range(3),{0,2},{},profile_cells=[0,1,2],halfspaces=h,curvature_limits_per_m={1:80.},on_step=lambda s,d:saved.append(s))
        self.assertTrue(r['qualifiedCourse']);self.assertGreaterEqual(len(saved),2)
        self.assertFalse(saved[0]['qualifiedCourse']);self.assertGreater(saved[0]['proof']['curvaturePhaseSlack_per_m'],0.)
        self.assertLessEqual(r['steps'][-1]['proof']['maximumActualGeometricCurvature_per_m'],80.+1e-4)
        self.assertTrue(np.array_equal(c,original));self.assertTrue(np.array_equal(delta[[0,2]],np.zeros((2,3))))

    def test_true_geometric_conflict_keeps_a_diagnostic_and_never_qualifies(self):
        c,h=self.fixture(.007)
        delta,r=fit_course(c,[[0,1],[1,2]],range(3),{0,2},{},profile_cells=[0,1,2],halfspaces=h,curvature_limits_per_m={1:80.},maximum_steps=3)
        self.assertIsNotNone(delta);self.assertFalse(r['qualifiedCourse']);self.assertEqual(len(r['steps']),3)
        self.assertAlmostEqual(delta[1,0],.007,places=7);self.assertTrue(np.array_equal(delta[[0,2]],np.zeros((2,3))))
        self.assertTrue(all(s['proof']['maximumHalfspaceViolation_m']<=1e-8 for s in r['steps']))

    def test_contradictory_hard_constraints_cannot_be_paid_off_by_curvature_slack(self):
        c,h=self.fixture(.007);h[1]['minimum_m']=-.006
        delta,r=fit_course(c,[[0,1],[1,2]],range(3),{0,2},{},profile_cells=[0,1,2],halfspaces=h,curvature_limits_per_m={1:80.})
        self.assertIsNone(delta);self.assertFalse(r['qualifiedCourse']);self.assertEqual(len(r['steps']),1)


if __name__=='__main__':unittest.main()
