import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.declared_source_fin import remove_declared_fin


class DeclaredSourceFinTests(unittest.TestCase):
    def mesh(self,*,aliases=False,repeat=False):
        p=np.array([[1.,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]])
        f=np.array([[0,2,4],[2,1,4],[1,3,4],[3,0,4],[2,0,5],[1,2,5],[3,1,5],[0,3,5],[0,2,1]])
        if repeat:f=np.vstack([f,f[-1]])
        if aliases:
            p=np.vstack([p,p[[0,2,1]]]);f[-1]=[6,7,8]
        return Mesh(p,f,'explicit-fixture')

    def test_internal_fin_removal_has_no_new_boundary(self):
        m=self.mesh();keep,proof=remove_declared_fin(m,8)
        np.testing.assert_array_equal(keep,np.arange(8));self.assertEqual(proof['newBoundaryEdges'],[])
        self.assertFalse(proof['anatomicalArtifactIdentityInferred'])

    def test_real_surface_face_cannot_be_deleted(self):
        with self.assertRaises(ValueError):remove_declared_fin(self.mesh(),0)

    def test_ambiguous_multiple_fins_rejected(self):
        with self.assertRaises(ValueError):remove_declared_fin(self.mesh(repeat=True),8)

    def test_exact_source_aliases_are_not_moved_or_welded(self):
        m=self.mesh(aliases=True);before=m.vertices.copy();keep,_=remove_declared_fin(m,8)
        np.testing.assert_array_equal(m.vertices,before);self.assertEqual(len(keep),8)
