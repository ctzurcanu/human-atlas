import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.algorithms.arterial_side_branch import make_side_branch

class SideBranchTests(unittest.TestCase):
    def test_real_portal_retains_parent_and_open_daughter(self):
        angles=np.arange(8)*np.pi/4;v=np.concatenate([np.c_[.001*np.cos(angles),.001*np.sin(angles),np.full(8,z)]for z in [0,.003,.006]])
        faces=[]
        for k in range(2):
            for j in range(8):
                a=k*8+j;b=k*8+(j+1)%8;faces.extend([[a,b,b+8],[a,b+8,a+8]])
        f=np.array(faces);normal=np.tile(np.c_[np.cos(angles),np.sin(angles),np.zeros(8)],(3,1));attrs=dict(normals=np.rint(normal*32767).astype(np.int16),normalSigns=np.ones(len(v),np.int8));m=Mesh(v,f,'test:artery',attrs)
        tri=v[f[0]];n=np.cross(tri[1]-tri[0],tri[2]-tri[0]);n/=np.linalg.norm(n)
        patched,branch,joined,proof=make_side_branch(m,0,tri.mean(0)+n*.006,n,radius_m=.0001)
        np.testing.assert_array_equal(patched.vertices[:len(v)],v)
        np.testing.assert_array_equal(patched.attributes['normals'][:len(v)],attrs['normals'])
        np.testing.assert_array_equal(patched.vertices[proof['parentPortalRawIds']],branch.vertices[proof['daughterPortalRawIds']])
        np.testing.assert_array_equal(patched.attributes['normals'][proof['parentPortalRawIds']],branch.attributes['normals'][proof['daughterPortalRawIds']])
        check=diagnostics(joined,self_intersections=True)
        self.assertEqual(check['self_intersection_pairs'],0);self.assertEqual(check['nonmanifold_edges'],0)
        e,c=np.unique(np.sort(joined.faces[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1),axis=0,return_counts=True)
        self.assertEqual(int((c==1).sum()),8+8+6)
        self.assertEqual(proof['lumenDisksCreated'],0);self.assertLess(proof['maximumRelativeSectionAreaError'],1e-10)

    def test_no_silent_calibre_reduction(self):
        m=Mesh(np.array([[0,0,0],[.0001,0,0],[0,.0001,0]]),np.array([[0,1,2]]),'test:artery',dict(normals=np.tile([0,0,32767],(3,1)).astype(np.int16),normalSigns=np.ones(3,np.int8)))
        with self.assertRaisesRegex(ValueError,'calibre'):make_side_branch(m,0,[0,0,.01],[0,0,1],radius_m=.00015)

if __name__=='__main__':unittest.main()
