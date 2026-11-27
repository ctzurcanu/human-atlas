import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_boundary_slit import zip_source_slit

def tube(gap=1e-5):
    ring=np.array([[.002,0.],[0.,.002],[-.002,0.],[0.,-.002]])
    v=np.array([[x,y,z] for z in np.arange(7)*.001 for x,y in ring])
    duplicate={8:28,12:29,16:30}
    v=np.r_[v,v[[8,12,16]]];v[[28,30],1]-=gap
    faces=[]
    for layer in range(6):
        for j in range(4):
            a=layer*4+j;b=layer*4+(j+1)%4;c=b+4;d=a+4
            tris=[[a,b,c],[a,c,d]]
            if j==3:tris=[[duplicate.get(i,i) for i in t] for t in tris]
            faces+=tris
    f=np.array(faces);tri=v[f];fn=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);normal=np.zeros_like(v)
    for j in range(3):np.add.at(normal,f[:,j],fn)
    normal/=np.linalg.norm(normal,axis=1)[:,None]
    attrs={'normals':np.rint(normal*32767).astype(np.int16),'normalSigns':np.ones(len(v),np.int8),
           'colors':np.arange(len(v)*3,dtype=np.uint8).reshape(-1,3)}
    return Mesh(v,f,'test:metres',attrs)

class SourceBoundarySlitTests(unittest.TestCase):
    def test_only_declared_shore_controls_move_and_real_mouths_stay_open(self):
        m=tube();q,p=zip_source_slit(m,[4,8,12,16,20],[4,28,29,30,20],maximum_motion_m=3e-5)
        np.testing.assert_array_equal(np.flatnonzero(np.any(q.vertices!=m.vertices,axis=1)),[28,30])
        np.testing.assert_array_equal(q.faces,m.faces)
        np.testing.assert_array_equal(q.attributes['colors'],m.attributes['colors'])
        self.assertEqual(p['genuineSourceMouthsRetained'],2)
        self.assertEqual(p['facesCreated'],0);self.assertEqual(p['capsCreated'],0)

    def test_wide_gap_is_not_a_source_seam(self):
        with self.assertRaisesRegex(ValueError,'narrow'):
            zip_source_slit(tube(.0003),[4,8,12,16,20],[4,28,29,30,20],maximum_motion_m=3e-5)

    def test_nonmatching_endpoint_rejected(self):
        with self.assertRaisesRegex(ValueError,'endpoints'):
            zip_source_slit(tube(),[4,8,12,16,20],[4,28,29,30,19],maximum_motion_m=3e-5)

    def test_incomplete_boundary_not_filled(self):
        with self.assertRaisesRegex(ValueError,'boundary|source crack'):
            zip_source_slit(tube(),[4,8,16,20],[4,28,30,20],maximum_motion_m=3e-5)

if __name__=='__main__':unittest.main()
