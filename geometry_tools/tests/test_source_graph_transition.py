import unittest
import numpy as np
from geometry_tools.algorithms.source_graph_transition import fixed_interface_weights

class SourceGraphTransitionTests(unittest.TestCase):
 def strip(self):
  v=np.array([[x,y,0]for x in (0.,.01,.02,.03)for y in (0.,.001)])
  f=np.array([[2*i,2*i+1,2*i+2]for i in range(3)]+[[2*i+1,2*i+3,2*i+2]for i in range(3)])
  return v,f,np.repeat(np.arange(4),2)
 def test_fixed_port_and_source_path_transition(self):
  v,f,c=self.strip();delta=np.tile([.005,0,0],(4,1));w,e=fixed_interface_weights(v,f,c,[0],delta)
  np.testing.assert_allclose(e['distance_m'],[0,.01,.02,.03]);np.testing.assert_allclose(w,[0,.5,1,1],atol=1e-12)
 def test_large_displacement_requires_longer_source_transition(self):
  v,f,c=self.strip();delta=np.tile([.01,0,0],(4,1));w,e=fixed_interface_weights(v,f,c,[0],delta)
  self.assertAlmostEqual(e['required_transition_m'][2],.04);self.assertAlmostEqual(w[2],.5)
 def test_nearby_disconnected_source_is_not_tethered(self):
  v,f,c=self.strip();v=np.r_[v,[[0,0,.00001],[.001,0,.00001],[0,.001,.00001]]];f=np.r_[f,[[8,9,10]]];c=np.r_[c,[4,4,4]]
  w,e=fixed_interface_weights(v,f,c,[0],np.zeros((5,3)));self.assertTrue(np.isinf(e['distance_m'][4]));self.assertEqual(w[4],1)
 def test_invalid_transition_is_rejected(self):
  v,f,c=self.strip()
  with self.assertRaises(ValueError):fixed_interface_weights(v,f,c,[0],np.zeros((4,3)),transition_m=0)

if __name__=='__main__':unittest.main()
