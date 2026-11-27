"""A flat flip transition joins two triangulations without adding tissue."""
import unittest
import numpy as np
from geometry_tools.algorithms.tetrahedral_orientation import orient_complex
from geometry_tools.implicit_cells import Arrangement
from geometry_tools.model import Mesh
from geometry_tools.audit import audit


class TetrahedralOrientationTests(unittest.TestCase):
    def fixture(self):
        p=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0],[.4,.4,1],[.4,.4,-1]])
        t=np.array([[0,1,2,5],[0,2,3,5],[0,1,3,4],[1,2,3,4],[0,1,2,3]])
        det=np.linalg.det(p[t[:4,1:]]-p[t[:4,:1]])
        for i in np.flatnonzero(det<0):t[i]=t[i,[0,2,1,3]]
        boundary=[]
        for tet in t[:4]:
            for row in ((1,2,3),(0,3,2),(0,1,3),(0,2,1)):
                f=tet[list(row)]
                if np.any(p[f,2]) and not (4 in f and 1 in f and 3 in f) and not (5 in f and 0 in f and 2 in f):
                    boundary.append(f)
        return p,t,np.array(boundary)

    def test_flat_transition_uses_adjacency_and_cancels_all_internal_faces(self):
        p,t,boundary=self.fixture()
        result,proof=orient_complex(p,t,boundary)
        self.assertEqual(proof['flat_interface_glue_tetrahedra'],[4])
        a=Arrangement(p,np.full((len(p),1),-1.),lambda x:'one')
        v,s,c=a.build(result)
        mesh=Mesh(v,next(iter(s.values())),'synthetic:metres')
        checked=audit(mesh,require_closed=True,expected_components=1,expected_genus=[0])
        self.assertTrue(checked['passed'],checked['failures'])
        self.assertAlmostEqual(checked['topology']['oriented_triangle_integral_m3'],2/3)
        reversed=t.copy();reversed[4]=reversed[4,[0,2,1,3]]
        corrected,_=orient_complex(p,reversed,boundary)
        np.testing.assert_array_equal(result,corrected)

    def test_real_inversion_is_not_excused_as_a_flat_transition(self):
        p,t,boundary=self.fixture();p[4,2]=-2
        with self.assertRaisesRegex(ValueError,'inverted'):
            orient_complex(p,t,boundary)


if __name__=='__main__':
    unittest.main()
