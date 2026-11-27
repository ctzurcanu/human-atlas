import unittest
import numpy as np
from geometry_tools.algorithms.source_port_candidates import source_port_candidates

class SourcePortCandidatesTests(unittest.TestCase):
 def port(self,name,x,moving=False,system='venous'):
  return dict(partId=name,system=system,moving=moving,source=np.array([[x,.001,0],[x,0,.001],[x,-.001,0],[x,0,-.001]]))
 def test_complete_inventory_retains_stationary_connection(self):
  ports=[self.port('parent',0),self.port('branch',.000001)]
  self.assertEqual(list(source_port_candidates(ports)),[(0,1)])
  self.assertEqual(list(source_port_candidates(ports,include_unmoved=False)),[])
 def test_inventory_does_not_depend_on_previous_candidate_motion(self):
  ports=[self.port('a',0),self.port('b',.000001),self.port('c',.1),self.port('d',.100001)]
  expected=[(0,1),(2,3)];self.assertEqual(list(source_port_candidates(ports)),expected)
  ports[0]['moving']=True;self.assertEqual(list(source_port_candidates(ports)),expected)
  self.assertEqual(list(source_port_candidates(ports,include_unmoved=False)),[(0,1)])
 def test_source_flow_and_identity_filter_before_contour_comparison(self):
  ports=[self.port('a',0),self.port('a',0),self.port('artery',0,system='arterial')]
  self.assertEqual(list(source_port_candidates(ports)),[])
 def test_same_mesh_record_requires_distinct_explicit_components(self):
  ports=[dict(self.port('a',0),component=0),dict(self.port('a',.000001),component=1),dict(self.port('a',.000002),component=0)]
  self.assertEqual(list(source_port_candidates(ports,include_same_part=True)),[(0,1),(1,2)])
 def test_same_part_without_component_ownership_is_rejected(self):
  with self.assertRaises(ValueError):list(source_port_candidates([self.port('a',0),self.port('a',0)],include_same_part=True))

if __name__=='__main__':unittest.main()
