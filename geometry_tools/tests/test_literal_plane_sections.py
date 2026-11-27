import unittest
import numpy as np
from geometry_tools.algorithms.literal_plane_sections import literal_plane_sections

class LiteralSectionsTest(unittest.TestCase):
 def test_closed_cube_cut_uses_actual_edges(self):
  v=np.array([[x,y,z]for x in [-1,1]for y in [-1,1]for z in [-1,1]],float);f=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
  r=literal_plane_sections(v,f,[0,0,.2],[0,0,1]);self.assertTrue(r['closed_local_sections']);self.assertEqual(len(r['loops']),1);self.assertTrue(np.allclose(r['loops'][0][:,2],.2));self.assertEqual(r['caps_created'],0)
 def test_open_sheet_is_not_filled(self):
  r=literal_plane_sections([[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]],[[0,1,2],[0,2,3]],[0,0,0],[1,0,0]);self.assertFalse(r['closed_local_sections']);self.assertEqual(len(r['loops']),0);self.assertEqual(len(r['open_components']),1)
 def test_separate_loops_retain_literal_source_face_ownership(self):
  v=np.array([[x,y,z]for x in [-1,1]for y in [-1,1]for z in [-1,1]],float);f=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
  r=literal_plane_sections(np.r_[v,v+[3,0,0]],np.r_[f,f+8],[0,0,.2],[0,0,1]);self.assertEqual(len(r['loops']),2);owners=r['loop_source_faces'];self.assertTrue(all(all(i<12 for i in row)or all(i>=12 for i in row)for row in owners));self.assertFalse(set(owners[0])&set(owners[1]));self.assertEqual(set(sum(owners,[])),set(r['section_source_faces']))
 def test_edge_weights_replay_normal_preserving_motion_exactly(self):
  v=np.array([[x,y,z]for x in [-1,1]for y in [-1,1]for z in [-1,1]],float);f=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]]);r=literal_plane_sections(v,f,[0,.2,0],[0,1,0]);q=v.copy();q[:,[0,2]]+=np.random.default_rng(2).normal(0,.1,(8,2));after=literal_plane_sections(q,f,[0,.2,0],[0,1,0]);predicted=np.sum(q[r['loop_raw_edges'][0]]*r['loop_edge_weights'][0][:,:,None],axis=1);actual=after['loops'][0];self.assertLess(np.min(np.linalg.norm(predicted[:,None]-actual[None],axis=2),axis=1).max(),1e-12)
 def test_close_disconnected_sheets_are_not_welded(self):
  v=np.array([[-1,-1,0],[1,-1,0],[1,1,0],[-1,1,0]],float);v=np.r_[v,v+[0,0,1e-10]];f=np.array([[0,1,2],[0,2,3],[4,5,6],[4,6,7]]);r=literal_plane_sections(v,f,[0,0,0],[1,0,0]);self.assertEqual(len(r['open_components']),2);self.assertEqual(len(r['loops']),0)
if __name__=='__main__':unittest.main()
