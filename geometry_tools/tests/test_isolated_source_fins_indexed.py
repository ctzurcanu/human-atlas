import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.isolated_source_fins_indexed import isolated_source_fins_indexed as isolated_source_fins

class SourceFinTests(unittest.TestCase):
    def plane(self):
        vertices=np.array([[i,j,0.] for j in range(4) for i in range(4)])
        faces=[]
        for j in range(3):
            for i in range(3):
                a=j*4+i;faces.extend([[a,a+1,a+5],[a,a+5,a+4]])
        return vertices,np.array(faces)

    def test_only_dangling_quad_is_removed(self):
        v,f=self.plane();v=np.r_[v,[[1.,.5,.5],[.5,.5,.5]]]
        # Existing opposed wall faces share diagonal (0,5).
        faces=np.r_[f,[[0,5,16],[0,16,17]]]
        m=Mesh(v,faces,'test');keep,p=isolated_source_fins(m)
        np.testing.assert_array_equal(keep,np.arange(len(f)))
        self.assertEqual(p['nonmanifoldEdgesBefore'],1);self.assertEqual(p['nonmanifoldEdgesAfter'],0)
        self.assertEqual(p['newBoundaryEdges'],0)

    def test_regular_open_quad_is_retained(self):
        v,f=self.plane();m=Mesh(v,f,'test');keep,p=isolated_source_fins(m)
        np.testing.assert_array_equal(keep,np.arange(len(f)));self.assertEqual(p['patches'],[])

    def test_three_face_branch_is_not_a_two_triangle_fin(self):
        v,f=self.plane();v=np.r_[v,[[1.,.5,.5],[.5,.5,.5],[.5,.8,.9]]]
        m=Mesh(v,np.r_[f,[[0,5,16],[0,16,17],[17,16,18]]],'test')
        keep,p=isolated_source_fins(m)
        self.assertEqual(len(keep),len(m.faces));self.assertEqual(p['patches'],[])

if __name__=='__main__':unittest.main()
