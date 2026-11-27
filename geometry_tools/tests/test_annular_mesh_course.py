import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.annular_mesh_course import annular_mesh_course


class AnnularCourseTests(unittest.TestCase):
    def mesh(self):
        n=12;angles=np.arange(n)*2*np.pi/n;ring=np.c_[.001*np.cos(angles),.001*np.sin(angles),np.zeros(n)]
        centres=np.array([[0,0,0],[.003,0,.01],[.008,0,.02]])
        v=np.concatenate([ring+c for c in centres]);f=[]
        for k in range(2):
            for i in range(n):
                a=k*n+i;b=k*n+(i+1)%n;c=b+n;d=a+n;f.extend([[a,b,c],[a,c,d]])
        return Mesh(v,np.array(f),'test:metres'),centres

    def test_recovers_the_existing_course_without_assuming_face_order(self):
        m,c=self.mesh();m=Mesh(m.vertices,m.faces[::-1],m.frame)
        course,proof=annular_mesh_course(m,[0,0,-1])
        np.testing.assert_allclose(course['centres'],c,atol=1e-15)
        np.testing.assert_allclose(course['radius_m'],.001,atol=1e-15)
        self.assertTrue(proof['everySourceFaceAccountedFor'])

    def test_a_missing_side_face_is_not_silently_turned_into_a_complete_guide(self):
        m,_=self.mesh()
        with self.assertRaises(ValueError):annular_mesh_course(Mesh(m.vertices,m.faces[1:],m.frame),[0,0,-1])


if __name__=='__main__':unittest.main()
