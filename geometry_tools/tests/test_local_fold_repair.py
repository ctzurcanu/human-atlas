"""Selective group repair must retain pre-existing defects as explicit evidence."""
import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.algorithms.local_fold_repair import repair_folds

def fixture():
 pair=np.array([[-1.,0,0],[1,0,0],[0,1,0],[0,.2,-1],[0,.2,1],[0,-.5,0]])*1e-5
 v=np.concatenate([pair,pair+[.001,0,0]])
 f=np.arange(len(v)).reshape(-1,3);t=v[f];n=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]);n/=np.linalg.norm(n,axis=1)[:,None]
 return Mesh(v,f,'fixture',{'normals':np.repeat(np.rint(n*32767).astype('<i2'),3,axis=0),'normalSigns':np.ones(len(v),dtype='i1')})

class FoldTests(unittest.TestCase):
 def test_selective_target_does_not_claim_all_defects_repaired(self):
  m=fixture();before=m.geometry_sha
  result,proof=repair_folds(m,[[0,1],[2,3]],target_pairs=[[0,1]],maximum_motion_m=.0001,wall_seconds=2)
  np.testing.assert_array_equal(result.faces,m.faces);np.testing.assert_array_equal(result.vertices[6:],m.vertices[6:])
  self.assertEqual(m.geometry_sha,before);self.assertEqual(proof['remainingTargetPairs'],0);self.assertFalse(proof['allSourceDefectsRemoved'])
  d=diagnostics(result,include_pairs=True);self.assertEqual(d['source_face_pairs'],[[2,3]])

 def test_target_and_owner_scope_are_explicit(self):
  m=fixture()
  with self.assertRaisesRegex(ValueError,'target'):
   repair_folds(m,[[0,1],[2,3]],target_pairs=[[0,2]])
  with self.assertRaisesRegex(ValueError,'owner'):
   repair_folds(m,[[0,1],[2,3]],face_owners=[0])
  with self.assertRaisesRegex(ValueError,'Unresolved'):
   repair_folds(m,[[0,1],[2,3]],locked=range(len(m.vertices)))

if __name__=='__main__':unittest.main()
