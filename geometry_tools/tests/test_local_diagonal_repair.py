import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.algorithms.local_diagonal_repair import repair_diagonals


class LocalDiagonalTests(unittest.TestCase):
 def fixture(self,depth=.0002):
  p=np.array([[0,0,0],[1,0,0],[1,1,depth],[0,1,0],[.45,.45,depth*.475],[.55,.45,depth*.475],[.5,.55,depth*.475]])
  return Mesh(p,np.array([[0,1,2],[0,2,3],[4,5,6]]),'analytical-metres',dict(uvs=p[:,:2]))
 def test_diagonal_removes_crossing_without_moving_or_capping(self):
  m=self.fixture();q,proof=repair_diagonals(m);self.assertGreater(proof['initialCrossings'],0);self.assertEqual(proof['finalCrossings'],0);np.testing.assert_array_equal(q.vertices,m.vertices);np.testing.assert_array_equal(q.attributes['uvs'],m.attributes['uvs']);self.assertEqual(len(q.faces),len(m.faces));self.assertEqual(topology(q)['open_edges'],topology(m)['open_edges']);np.testing.assert_array_equal(m.faces,[[0,1,2],[0,2,3],[4,5,6]])
 def test_thick_patch_cannot_be_silently_changed(self):
  with self.assertRaises(ValueError):repair_diagonals(self.fixture(depth=.02))
 def test_explicit_partial_repair_retains_all_unresolved_crossings_without_new_pairs(self):
  first,second=self.fixture(),self.fixture(depth=.02)
  m=Mesh(np.r_[first.vertices,second.vertices+[3,0,0]],np.r_[first.faces,second.faces+len(first.vertices)],first.frame)
  with self.assertRaises(ValueError):repair_diagonals(m)
  q,proof=repair_diagonals(m,allow_partial=True)
  self.assertTrue(proof['complete']);self.assertFalse(proof['allIntrinsicCrossingsResolved'])
  self.assertGreater(proof['initialCrossings'],proof['finalCrossings']);self.assertGreater(proof['finalCrossings'],0)
  self.assertEqual(proof['newIntrinsicCrossingPairs'],0)
  from geometry_tools.open3d_backend import diagnostics
  original=diagnostics(m,include_pairs=True);after=diagnostics(q,include_pairs=True)
  self.assertLess(set(map(tuple,after['source_face_pairs'])),set(map(tuple,original['source_face_pairs'])))
  np.testing.assert_array_equal(m.vertices,q.vertices)
  np.testing.assert_array_equal(m.faces[3:],q.faces[3:])
 def test_timeout_never_becomes_partial_acceptance(self):
  with self.assertRaises(TimeoutError):repair_diagonals(self.fixture(),wall_seconds=0,allow_partial=True)
 def test_protected_contact_face_cannot_change(self):
  with self.assertRaises(ValueError):repair_diagonals(self.fixture(),allow_partial=True,protected_face_ids=[0])
  with self.assertRaises(ValueError):repair_diagonals(self.fixture(),protected_face_ids=[3])

if __name__=='__main__':unittest.main()
