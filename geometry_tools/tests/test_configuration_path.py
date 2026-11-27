import unittest
import numpy as np
from geometry_tools.algorithms.configuration_path import configuration_path,curved_primary_configuration_path


class ConfigurationPathTest(unittest.TestCase):
 def test_smooth_route_can_prefer_a_more_expensive_station(self):
  states=[np.array([[[0,0]]]),np.array([[[0,0]],[[.01,0]]]),np.array([[[0,0]]])]
  path,proof=configuration_path(states,[[0],[.0001,0],[0]],[[(0,0,1)],[(0,0,1)]])
  self.assertEqual(path,[0,0,0]);self.assertTrue(proof['feasible'])
 def test_no_continuous_route_is_explicit(self):
  path,proof=configuration_path([np.array([[[0,0]]]),np.array([[[.1,0]]])],[[0],[0]],[[(0,0,1)]])
  self.assertIsNone(path);self.assertFalse(proof['feasible'])
 def test_minimax_diagnostic_measures_required_curvature_without_waiving_the_bound(self):
  states=[np.array([[[0.,0.]]]),np.array([[[.1,0.]],[[.2,0.]]]),np.array([[[0.,0.]]])]
  path,proof=curved_primary_configuration_path(states,[[0.],[100.,0.],[0.]],[[(0,0,1.)],[(0,0,1.)]],[0.,1.,2.],maximum_primary_step_m=1.,maximum_curvature_per_m=.01,measure_minimum_maximum_curvature=True)
  self.assertEqual(path,[0,0,0]);self.assertTrue(proof['diagnosticOnly']);self.assertAlmostEqual(proof['maximumPrimaryCurvature_per_m'],.2/1.01)
 def test_second_order_search_rejects_cheaper_zigzag(self):
  states=[np.array([[[0.,0.]]]),np.array([[[0.,0.]],[[.1,0.]]]),np.array([[[0.,0.]]])]
  path,proof=curved_primary_configuration_path(states,[[0.],[1.,0.],[0.]],[[(0,0,1.)],[(0,0,1.)]],[0.,1.,2.],maximum_primary_step_m=1.,maximum_curvature_per_m=.1)
  self.assertEqual(path,[0,0,0]);self.assertEqual(proof['maximumPrimaryCurvature_per_m'],0.)
 def test_second_order_impossible_bend_is_not_claimed_smooth(self):
  states=[np.array([[[0.,0.]]]),np.array([[[1.,0.]]]),np.array([[[0.,0.]]])]
  path,proof=curved_primary_configuration_path(states,[[0.],[0.],[0.]],[[(0,0,1.)],[(0,0,1.)]],[0.,.01,.02],maximum_primary_step_m=2.,maximum_curvature_per_m=1.)
  self.assertIsNone(path);self.assertFalse(proof['feasible'])
