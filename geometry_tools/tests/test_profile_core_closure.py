import unittest
import numpy as np
from geometry_tools.algorithms.profile_core_transport import incident_profile_core_closure


class ProfileCoreClosureTest(unittest.TestCase):
 def test_complete_connected_core_stops_at_other_sections(self):
  # Sections 0 and 1 bound a long branch core. Section 2 and its separate
  # branch are unrelated. Coordinates never enter the connectivity query.
  profiles=[np.array([0,1,2]),np.array([3,4,5]),np.array([6,7,8])]
  faces=np.array([[0,1,9],[1,9,10],[9,10,11],[10,11,3],[3,4,12],[6,7,13]])
  r=incident_profile_core_closure(14,faces,profiles,[0])
  self.assertEqual(np.flatnonzero(r['core_mask']).tolist(),[9,10,11])
  self.assertEqual(r['components'][0]['boundaryProfileIds'],[0,1])
  self.assertFalse(r['near_points_welded'])

 def test_empty_selection_and_multiple_incident_patches(self):
  profiles=[[0,1,2],[3,4,5]];faces=[[0,1,6],[0,2,7],[3,4,8]]
  self.assertFalse(incident_profile_core_closure(9,faces,profiles,[])['core_mask'].any())
  self.assertEqual(np.flatnonzero(incident_profile_core_closure(9,faces,profiles,[0])['core_mask']).tolist(),[6,7])

 def test_rejects_overlapping_sections_and_unknown_controls(self):
  with self.assertRaises(ValueError):incident_profile_core_closure(6,[[0,1,2]],[[0,1,2],[2,3,4]],[0])
  with self.assertRaises(ValueError):incident_profile_core_closure(6,[[0,1,6]],[[0,1,2]],[0])
  with self.assertRaises(ValueError):incident_profile_core_closure(6,[[0,1,2]],[[0,1,2]],[1])


if __name__=='__main__':unittest.main()
