import unittest
import numpy as np
from geometry_tools.algorithms.annular_member_pose import AnnularMemberPose,transport_rigid_cells


class AnnularMemberPoseTest(unittest.TestCase):
 def test_proper_pose_tapers_to_literal_source(self):
  s=dict(y_m=[.19,.20,.21],ySupport_m=[.17,.19,.21,.24],sourceArteryCentreXZ_m=[[0,0]]*3,targetXZRoll=[[.003,0,.2],[.004,0,.3],[.003,0,.2]],maximumChange_m=.008)
  field=AnnularMemberPose(s)
  for y in [.16,.17,.18,.20,.23,.24,.25]:
   R,t=field.pose(y,np.eye(3),np.zeros(3));self.assertAlmostEqual(np.linalg.det(R),1);self.assertTrue(np.allclose(R.T@R,np.eye(3)))
   if y<=.17 or y>=.24:self.assertTrue(np.array_equal(R,np.eye(3)));self.assertTrue(np.array_equal(t,np.zeros(3)))
  R,t=field.pose(.20,np.eye(3),np.zeros(3));self.assertTrue(np.allclose(t,[.004,0,0]))
 def test_literal_cells_and_fixed_mouth_remain_rigid(self):
  v=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[0,.5,0],[1,.5,0]],float)*.001
  faces=np.array([[0,1,5],[0,5,4],[4,5,3],[4,3,2]]);cells=np.array([0,0,1,1,2,3]);R=np.array([[0,-1,0],[1,0,0],[0,0,1]])
  q,e=transport_rigid_cells(v,faces,cells,{1:(R,np.array([.001,0,0]))},{0})
  self.assertTrue(np.array_equal(q[:2],v[:2]));self.assertTrue(np.allclose(q[2:4],v[2:4]@R.T+[.001,0,0]));self.assertLess(e['linearResidual'],1e-8)
  with self.assertRaises(ValueError):transport_rigid_cells(v,faces,cells,{0:(np.eye(3),np.array([.001,0,0]))},{0})
 def test_declared_source_component_end_is_an_exact_pose_anchor(self):
  field=AnnularMemberPose(dict(y_m=[.19,.20,.21],ySupport_m=[.17,.19,.21,.24],sourceArteryCentreXZ_m=[[0,0]]*3,targetXZRoll=[[.003,0,.2]]*3,maximumChange_m=.008,anchorYIntervals_m=[[.217,.219]],anchorTransition_m=.012))
  R,t=field.pose(.218,np.eye(3),np.zeros(3));self.assertTrue(np.array_equal(R,np.eye(3)));self.assertTrue(np.array_equal(t,np.zeros(3)))
  R,t=field.pose(.20,np.eye(3),np.zeros(3));self.assertGreater(np.linalg.norm(t),.002)


if __name__=='__main__':unittest.main()
