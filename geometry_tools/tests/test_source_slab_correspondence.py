import unittest
import numpy as np
from geometry_tools.algorithms.source_slab_correspondence import source_slab_correspondence
from geometry_tools.algorithms.literal_plane_sections import literal_plane_sections


def tubes_joined_only_above_slab():
 ring=np.array([[-.2,-.2],[.2,-.2],[.2,.2],[-.2,.2]])
 points=[];faces=[]
 for x in [-1.,1.]:
  base=len(points);points.extend(np.concatenate([np.c_[ring[:,0]+x,np.full(4,y),ring[:,1]]for y in [0,1,2,3]]))
  for row in range(3):
   for k in range(4):a=base+row*4+k;b=base+row*4+(k+1)%4;faces.extend([[a,b,b+4],[a,b+4,a+4]])
 faces.extend([[12,13,29],[12,29,28]])
 return np.array(points),np.array(faces)


class SourceSlabCorrespondenceTest(unittest.TestCase):
 def test_remote_source_join_does_not_mix_the_two_strands(self):
  v,f=tubes_joined_only_above_slab();a=literal_plane_sections(v,f,[0,1.2,0],[0,1,0]);b=literal_plane_sections(v,f,[0,1.8,0],[0,1,0]);self.assertEqual(len(a['loops']),2);r=source_slab_correspondence(v,f,1.2,1.8,a['loop_source_faces'],b['loop_source_faces']);self.assertEqual(len(r['links']),2);self.assertTrue(r['connectivityOutsideSlabExcluded']);self.assertFalse(any(c['explicitForkOrMerge']for c in r['components']))
 def test_source_profile_order_is_not_an_identity(self):
  v,f=tubes_joined_only_above_slab();a=literal_plane_sections(v,f,[0,1.2,0],[0,1,0]);b=literal_plane_sections(v,f,[0,1.8,0],[0,1,0]);r=source_slab_correspondence(v,f,1.2,1.8,a['loop_source_faces'],b['loop_source_faces'][::-1]);self.assertEqual({tuple(x[:2])for x in r['links']},{(0,1),(1,0)})
