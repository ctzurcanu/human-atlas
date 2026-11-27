import unittest
import numpy as np
from geometry_tools.algorithms.source_open_fin import remove_source_open_fin
from geometry_tools.model import Mesh

def fixture(radius=2):
    points=np.asarray([[x,y,0.] for x in range(-radius,radius+1) for y in range(-radius,radius+1)],float)
    lookup={tuple(p[:2]):i for i,p in enumerate(points)};faces=[]
    for x in range(-radius,radius):
        for y in range(-radius,radius):
            a,b,c,d=[lookup[p] for p in [(x,y),(x+1,y),(x+1,y+1),(x,y+1)]]
            faces.extend([[a,b,c],[a,c,d]])
    fin=[lookup[(1,0)],lookup[(0,0)],lookup[(0,1)]]
    return Mesh(points,np.asarray(faces+[fin]),'test:metres')

class SourceFinTests(unittest.TestCase):
    def test_extra_triangle_removed_and_wall_literal(self):
        m=fixture();before=m.vertices.copy();keep,proof=remove_source_open_fin(m,len(m.faces)-1)
        np.testing.assert_array_equal(keep,np.arange(len(m.faces)-1));np.testing.assert_array_equal(before,m.vertices)
        self.assertEqual(proof['newFaces'],0);self.assertTrue(proof['allOtherSourceBoundaryEdgesLiteral'])
    def test_true_wall_face_rejected(self):
        with self.assertRaises(ValueError):remove_source_open_fin(fixture(),0)
    def test_external_boundary_attachment_rejected(self):
        m=fixture(1)
        with self.assertRaisesRegex(ValueError,'external source-mouth'):remove_source_open_fin(m,len(m.faces)-1)
    def test_inconsistent_retained_winding_rejected(self):
        m=fixture();a,b=m.faces[-1,:2]
        ids=[i for i,f in enumerate(m.faces[:-1]) if a in f and b in f];self.assertEqual(len(ids),2)
        f=m.faces.copy();f[ids[0]]=f[ids[0],[0,2,1]];bad=Mesh(m.vertices,f,m.frame)
        with self.assertRaisesRegex(ValueError,'opposed'):remove_source_open_fin(bad,len(f)-1)

if __name__=='__main__':unittest.main()
