import unittest
import numpy as np
from geometry_tools.algorithms.section_profile_clearance import AnteriorSectionProfileClearance
from geometry_tools.tests.test_source_annular_refinement import tube

class SectionProfileClearanceTests(unittest.TestCase):
 def test_intruding_profile_moves_rigidly_to_anterior_side(self):
  v,f=tube([0,.02],radius=.001);oracle=AnteriorSectionProfileClearance(v,f);theta=np.arange(6)*np.pi/3;p=np.c_[.0004*np.cos(theta),np.full(6,.01),.0012+.0004*np.sin(theta)];shift,evidence=oracle.propose(p);self.assertGreater(shift[2],0);self.assertEqual(shift[0],0);self.assertEqual(shift[1],0);self.assertGreaterEqual(evidence['proposedMinimumGap_m'],.00005);np.testing.assert_allclose(np.linalg.norm(p[:,None]-p[None,:],axis=2),np.linalg.norm((p+shift)[:,None]-(p+shift)[None,:],axis=2),atol=1e-16)
 def test_clear_remote_profile_stays_fixed(self):
  v,f=tube([0,.02],radius=.001);oracle=AnteriorSectionProfileClearance(v,f);theta=np.arange(6)*np.pi/3;p=np.c_[.003+.0001*np.cos(theta),np.full(6,.01),.0001*np.sin(theta)];shift,_=oracle.propose(p);np.testing.assert_array_equal(shift,0)

if __name__=='__main__':unittest.main()
