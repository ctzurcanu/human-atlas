import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.model import Mesh
from geometry_tools.algorithms.anchored_section_pose import pose_anchored_sections
from geometry_tools.algorithms.source_profile_axis import source_profile_graph


class AnchoredSectionPoseTest(unittest.TestCase):
 def setUp(self):
  angle=np.arange(6)*2*np.pi/6
  self.p=np.array([[.002*np.cos(t),y,.002*np.sin(t)]for y in np.linspace(0,.06,7)for t in angle])
  self.f=np.array([[6*k+i,6*(k+1)+i,6*(k+1)+(i+1)%6]if z==0 else[6*k+i,6*(k+1)+(i+1)%6,6*k+(i+1)%6]for k in range(6)for i in range(6)for z in range(2)])
  self.mesh=Mesh(self.p,self.f,'test:metres')
  g=source_profile_graph(self.p,self.f,triangle_pairing=np.arange(len(self.f)).reshape(-1,2));self.ports=g['boundary_raw'];self.fixed=int(np.argmin([self.p[r,1].mean()for r in self.ports]));self.moving=1-self.fixed
  self.matrix=np.eye(4);self.matrix[:3,:3]=1.1*Rotation.from_rotvec([0,0,.15]).as_matrix();self.matrix[:3,3]=[.008,.008,.003]

 def test_fixed_and_moving_mouths_and_source_sections_survive(self):
  q,proof=pose_anchored_sections(self.mesh,self.matrix,fixed_ports=[self.fixed],moving_ports=[self.moving],unanchored_components='fixed')
  np.testing.assert_array_equal(q.faces,self.f);np.testing.assert_array_equal(q.vertices[self.ports[self.fixed]],self.p[self.ports[self.fixed]])
  ids=self.ports[self.moving];np.testing.assert_array_equal(q.vertices[ids],self.p[ids]@self.matrix[:3,:3].T+self.matrix[:3,3])
  self.assertTrue(all(x['maximumPairDistanceError_m']<1e-12 for x in proof['profiles']))
  self.assertAlmostEqual(proof['minimumSectionScale'],1.);self.assertAlmostEqual(proof['maximumSectionScale'],1.1)
  np.testing.assert_array_equal(self.mesh.vertices,self.p);self.assertTrue(proof['discreteAndNeighborAcceptancePending'])

 def test_explicit_complete_and_compatible_interface_policy_required(self):
  for fixed,moving in [([],[]),([0],[0]),([0],[3])]:
   with self.assertRaises(ValueError):pose_anchored_sections(self.mesh,self.matrix,fixed_ports=fixed,moving_ports=moving,unanchored_components='fixed')
  bad=self.matrix.copy();bad[:3,0]*=2
  with self.assertRaises(ValueError):pose_anchored_sections(self.mesh,bad,fixed_ports=[0],moving_ports=[1],unanchored_components='fixed')

if __name__=='__main__':unittest.main()
