import unittest
import numpy as np
from geometry_tools.algorithms.bounded_cell_fairing import fair_rigid_cells


class BoundedCellFairingTest(unittest.TestCase):
 def test_bounded_smooth_profile_chain_with_fixed_interfaces(self):
  c=np.c_[np.zeros(21),np.arange(21)*.002,np.zeros(21)];edges=np.c_[np.arange(20),np.arange(1,21)];d,e=fair_rigid_cells(c,edges,range(21),{0,1,19,20},{10:dict(translation_m=[.006,0])},profile_cells=range(21),curvature_per_m=60.)
  self.assertIsNotNone(d);self.assertTrue(e['feasible']);self.assertTrue(np.array_equal(d[[0,1,19,20]],np.zeros((4,3))));self.assertTrue(np.array_equal(d[:,1],np.zeros(21)));self.assertLessEqual(np.linalg.norm(d,axis=1).max(),.008+1e-8);self.assertLessEqual(e['maximumCurvatureBoundViolation_per_m'],1e-4);self.assertGreater(d[10,0],.002)
 def test_coupled_plane_moves_with_both_cells(self):
  c=np.array([[0,0,0],[.003,0,0],[0,.01,0],[.003,.01,0]]);edges=np.array([[0,2],[1,3]]);halfspaces=[dict(coefficients={0:[-1,0]},plane=(0,1),minimum_m=-.001),dict(coefficients={1:[1,0]},plane=(0,-1),minimum_m=-.001)];d,e=fair_rigid_cells(c,edges,range(4),{2,3},{0:dict(translation_m=[.005,0]),1:dict(translation_m=[.005,0])},halfspaces=halfspaces,plane_count=1)
  self.assertIsNotNone(d);self.assertLessEqual(e['maximumHalfspaceViolation_m'],1e-8);self.assertGreater(.003+d[1,0]-d[0,0],.0009)
 def test_infeasible_protected_constraint_is_not_accepted(self):
  c=np.array([[0,0,0],[0,.01,0]]);d,e=fair_rigid_cells(c,[[0,1]],range(2),{0},{},halfspaces=[dict(coefficients={0:[1,0]},minimum_m=.001)])
  self.assertIsNone(d);self.assertFalse(e['feasible'])
 def test_conflicting_free_cell_planes_have_replayable_dual_certificate(self):
  captured=[];c=np.array([[0,0,0],[0,.01,0]])
  d,e=fair_rigid_cells(c,[[0,1]],range(2),{1},{},halfspaces=[dict(coefficients={0:[1,0]},minimum_m=.002),dict(coefficients={0:[-1,0]},minimum_m=.001)],infeasibility_callback=lambda *args:captured.append(args))
  self.assertIsNone(d);self.assertTrue(e['infeasibilityCertificate']['verifiedNumerically']);self.assertEqual(len(captured),1)
  A,b,z,free,curvature=captured[0]
  self.assertAlmostEqual(float(b@z),-1.,places=7);self.assertLess(np.linalg.norm(A.T@z,np.inf),1e-6)
  self.assertTrue(np.array_equal(free,[0]));self.assertEqual(len(curvature),0);self.assertGreater(float(z[0]),0.);self.assertGreater(float(z[1]),0.)
 def test_refinement_edges_do_not_silently_remove_original_curve_obligations(self):
  c=np.array([[0,0,0],[0,.01,0],[0,.02,0],[.001,.01,0]])
  d,e=fair_rigid_cells(c,[[0,1],[1,2],[1,3]],range(4),{0,2},{1:dict(translation_m=[.002,0])},profile_cells=[0,1,2],curvature_edges=[[0,1],[1,2]])
  self.assertIsNotNone(d);self.assertTrue(e['feasible']);self.assertEqual(e['curvatureConstraints'],1)
 def test_current_arc_check_does_not_confuse_reference_parameter_with_geometric_curvature(self):
  c=np.array([[0.,0.,0.],[0.,.01,0.],[0.,.02,0.]]);reference=c.copy();reference[1,0]=.004
  d,e=fair_rigid_cells(c,[[0,1],[1,2]],range(3),{0,2},{1:dict(translation_m=[.007,0.],weight=100)},profile_cells=[0,1,2],curvature_reference_centres=reference,curvature_limits_per_m={1:80.})
  self.assertIsNotNone(d);self.assertTrue(e['feasible']);self.assertTrue(e['geometricCurvatureQualified'])
  self.assertGreater(d[1,0],.004);self.assertGreater(2*d[1,0]/.01**2,80.)
  actual=2*d[1,0]/(.01**2+d[1,0]**2)
  self.assertAlmostEqual(e['maximumActualGeometricCurvature_per_m'],actual,places=5);self.assertLessEqual(actual,80.+1e-4)
 def test_elastic_curvature_is_diagnostic_and_does_not_relax_fixed_or_clearance_constraints(self):
  c=np.array([[0.,0.,0.],[0.,.01,0.],[0.,.02,0.]])
  halfspaces=[dict(coefficients={1:[1.,0.]},minimum_m=.007),dict(coefficients={1:[-1.,0.]},minimum_m=-.007)]
  d,e=fair_rigid_cells(c,[[0,1],[1,2]],range(3),{0,2},{},profile_cells=[0,1,2],halfspaces=halfspaces,curvature_reference_centres=c,curvature_limits_per_m={1:80.},elastic_curvature_phase=True)
  self.assertIsNotNone(d);self.assertTrue(e['diagnosticProposalFeasible']);self.assertFalse(e['feasible']);self.assertFalse(e['geometricCurvatureQualified'])
  self.assertGreater(e['curvaturePhaseSlack_per_m'],1.);self.assertTrue(np.array_equal(d[[0,2]],np.zeros((2,3))));self.assertAlmostEqual(d[1,0],.007,places=8);self.assertLess(e['maximumHalfspaceViolation_m'],1e-8)
 def test_explicit_larger_placement_domain_still_enforces_its_limit(self):
  c=np.c_[np.zeros(9),np.arange(9)*.02,np.zeros(9)];edges=np.c_[np.arange(8),np.arange(1,9)]
  d,e=fair_rigid_cells(c,edges,range(9),{0,8},{4:dict(translation_m=[.020,0],weight=100)},maximum_motion_m=.025,profile_cells=range(9))
  self.assertIsNotNone(d);self.assertTrue(e['feasible']);self.assertGreater(d[4,0],.015);self.assertLessEqual(np.linalg.norm(d,axis=1).max(),.025+1e-8)
  with self.assertRaises(ValueError):fair_rigid_cells(c,edges,range(9),{0,8},{4:dict(translation_m=[.030,0])},maximum_motion_m=.025)


if __name__=='__main__':unittest.main()
