import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.rooted_source_sleeve import verify_sleeve_correspondence

class SourceSleeve(unittest.TestCase):
    def fixture(self,delta=1e-7):
        angle=np.arange(6)*2*np.pi/6;root=np.c_[np.cos(angle),np.sin(angle),np.zeros(6)];tip=root+[0,0,1];copy=tip+[delta,0,0]
        vertices=np.r_[root,tip,copy];faces=[]
        for base in [6,12]:
            for i in range(6):
                j=(i+1)%6;faces.extend([[i,j,base+j],[i,base+j,base+i]])
        m=Mesh(vertices,np.array(faces),'test');words=np.rint(vertices/1e-6).astype(int)
        return m,words

    def test_all_oriented_faces_have_a_bounded_correspondence(self):
        m,words=self.fixture();p=verify_sleeve_correspondence(m,np.arange(12,24),np.arange(12),np.arange(6),words)
        self.assertTrue(p['orientedSourceFaceBijection']);self.assertLessEqual(p['symmetricSurfaceDistanceUpperBound_m'],1.01e-7)
        self.assertFalse(p['sourcePositionsModified'])

    def test_geometry_and_source_code_limits_are_independent(self):
        m,words=self.fixture(delta=2e-6)
        with self.assertRaises(ValueError):verify_sleeve_correspondence(m,np.arange(12,24),np.arange(12),np.arange(6),words)
        m,words=self.fixture();words[12,0]+=5
        with self.assertRaises(ValueError):verify_sleeve_correspondence(m,np.arange(12,24),np.arange(12),np.arange(6),words)

    def test_opposite_winding_is_not_a_duplicate(self):
        m,words=self.fixture();faces=m.faces.copy();faces[12:]=faces[12:,::-1];wrong=Mesh(m.vertices,faces,m.frame)
        with self.assertRaises(ValueError):verify_sleeve_correspondence(wrong,np.arange(12,24),np.arange(12),np.arange(6),words)

if __name__=='__main__':unittest.main()
