from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ta98_triangle_intersections import crossing,proven_non_crossing_pairs,sweep_hits

class TrianglePairProofTests(unittest.TestCase):
 def test_rejections_never_hide_actual_crossings(self):
  rng=np.random.default_rng(713);a=np.array([[0.,0.,0.],[.002,0.,0.],[0.,.002,0.]])
  b=rng.normal(0,.001,(3000,3,3));b[:1000,:,2]=0;b[1000:1500,0]=a[0];b[1000:1500,1]=a[1]
  for i,width in enumerate(np.geomspace(1e-13,1e-5,100)):
   b[1500+i]=np.array([[0,0,0],[.002,0,0],[.001,width,0]])
  n=np.cross(b[:,1]-b[:,0],b[:,2]-b[:,0]);n/=np.linalg.norm(n,axis=1)[:,None];proof=proven_non_crossing_pairs(a,b,np.array([0.,0.,1.]),n)
  actual=np.array([crossing(a,t)for t in b]);self.assertFalse(np.any(actual&proof));self.assertGreater(proof.sum(),500);self.assertTrue(actual.any())
 def test_overlapping_coplanar_faces_are_still_reported(self):
  a=np.array([[0.,0.,0.],[.002,0.,0.],[0.,.002,0.]])
  b=a.copy();b[2]=[.0005,.001,0];hit=sweep_hits(np.array([a]),np.array([b]));self.assertEqual(hit['properOrCoplanarIntersections'],1)
if __name__=='__main__':unittest.main()
