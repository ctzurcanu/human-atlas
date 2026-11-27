import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
import numpy as np
from geometry_tools.algorithms.adaptive_course_transport import refine_course
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

class CourseTransportTest(unittest.TestCase):
    def source(self):
        v=np.array([[-1,-1,-1],[1,-1,-1],[1,1,-1],[-1,1,-1],[-1,-1,1],[1,-1,1],[1,1,1],[-1,1,1]],float)*.01
        f=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],[1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        # Every corner has a separate attribute ID. Geometry is still closed.
        return v[f].reshape(-1,3),np.arange(36).reshape(-1,3)
    def test_identity_is_literal(self):
        v,f=self.source();q,g,n,p,e=refine_course(v,f,np.tile([1.,0.,0.],(len(v),1)),lambda x:x.copy())
        self.assertTrue(np.array_equal(q,v));self.assertTrue(np.array_equal(g,f));self.assertTrue(np.array_equal(p,np.arange(len(f))))
    def test_curved_map_conforms_across_split_normal_edges(self):
        v,f=self.source()
        def field(x):
            q=x.copy();q[:,2]+=30*x[:,0]**2;q[:,1]+=20*x[:,0]**2;return q
        q,g,n,p,e=refine_course(v,f,np.tile([1.,0.,0.],(len(v),1)),field,deviation_m=2e-5,max_faces=10000,max_rounds=12)
        self.assertTrue(np.array_equal(q[:len(v)],v));self.assertEqual(len(g),len(p));self.assertTrue(set(p)<=set(range(len(f))))
        mapped=field(q);result=topology(Mesh(mapped,g,'test:metres'));self.assertTrue(result['closed_oriented_manifold'],result)
        edges=np.unique(np.sort(np.concatenate([g[:,[0,1]],g[:,[1,2]],g[:,[2,0]]]),axis=1),axis=0)
        self.assertLessEqual(np.linalg.norm(field(q[edges].mean(1))-mapped[edges].mean(1),axis=1).max(),2e-5)
        self.assertLessEqual(np.linalg.norm(field(q[g].mean(1))-mapped[g].mean(1),axis=1).max(),2e-5)

if __name__=='__main__':unittest.main()
