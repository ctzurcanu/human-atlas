import itertools
import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.internal_opposed_facets import internal_opposed_facets

def tetra(vertices):
    v=np.asarray(vertices,float);faces=[]
    for ids in itertools.combinations(range(4),3):
        f=list(ids);other=next(i for i in range(4) if i not in f);t=v[f]
        if np.cross(t[1]-t[0],t[2]-t[0])@(v[other]-t[0])>0:f=[f[0],f[2],f[1]]
        faces.append(f)
    return v,np.asarray(faces)

class InternalFacetTests(unittest.TestCase):
    def test_shared_internal_tetra_face_is_removed_with_solid_preserved(self):
        a,f=tetra([[0,0,0],[0,1,0],[0,0,1],[-1,0,0]])
        b,g=tetra([[0,0,0],[0,1,0],[0,0,1],[1,0,0]])
        m=Mesh(np.r_[a,b],np.r_[f,g+4],'fixture:metres')
        keep,proof=internal_opposed_facets(m)
        self.assertEqual(len(keep),6)
        self.assertEqual(len(proof['removedSourceFaceIds']),2)
        self.assertTrue(proof['retainedEnvelopeClosed'])
        self.assertLess(proof['maximumSourceWindingDifference'],1e-10)
        np.testing.assert_array_equal(m.faces,np.r_[f,g+4])

    def test_external_opposed_sheet_cannot_be_called_an_internal_partition(self):
        v,f=tetra([[0,0,0],[0,1,0],[0,0,1],[1,0,0]])
        points=np.r_[v,[[3,0,0],[3,1,0],[3,0,1]]]
        m=Mesh(points,np.r_[f,[[4,5,6],[6,5,4]]],'fixture:metres')
        with self.assertRaisesRegex(ValueError,'not inside'):
            internal_opposed_facets(m)

    def test_a_normal_closed_surface_without_opposed_faces_is_retained(self):
        v,f=tetra([[0,0,0],[0,1,0],[0,0,1],[1,0,0]])
        with self.assertRaisesRegex(ValueError,'No exactly opposed'):
            internal_opposed_facets(Mesh(v,f,'fixture:metres'))

if __name__=='__main__':unittest.main()
