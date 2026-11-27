import unittest
import numpy as np
from geometry_tools.algorithms.source_annulus_wall import reconstruct_declared_sidewall
from geometry_tools.algorithms.partial_annulus_closure import complete_five_quad_collar
from geometry_tools.algorithms.tube_continuity import continuity_audit


class SourceAnnulusWallTests(unittest.TestCase):
    def source(self):
        t=np.arange(6)*np.pi/3
        p=np.r_[np.c_[np.cos(t),np.sin(t),np.zeros(6)],np.c_[np.cos(t),np.sin(t),np.ones(6)]]*.001
        faces=[]
        for i in range(5):faces.extend([[i,(i+1)%6,i+6],[(i+1)%6,(i+1)%6+6,i+6]])
        return p,np.array(faces)

    def test_lateral_strip_has_two_open_ports_and_literal_controls(self):
        p,f=self.source();before=p.copy();added,proof=reconstruct_declared_sidewall(p,f,range(10),range(6))
        np.testing.assert_array_equal(p,before)
        self.assertEqual(added.shape,(2,3))
        self.assertTrue(proof['emittedGapFacesPermitted'])
        self.assertFalse(proof['computationOnly'])
        self.assertFalse(proof['rootOrEndCap'])
        self.assertTrue(continuity_audit(p,np.r_[f,added],expected_components=1,expected_mouths=2,expected_euler=0)['topological_continuity_pass'])
        for tri in added:self.assertEqual(len(np.unique(p[tri,2])),2)

    def test_computational_and_anatomical_permissions_are_distinct(self):
        p,f=self.source();added,proof=complete_five_quad_collar(p,f,range(10),range(6))
        self.assertTrue(proof['computationOnly'])
        self.assertFalse(proof['emittedGapFacesPermitted'])
        actual,_=reconstruct_declared_sidewall(p,f,range(10),range(6))
        np.testing.assert_array_equal(actual,added)

    def test_unowned_shape_or_wrong_root_is_rejected(self):
        p,f=self.source()
        with self.assertRaises(ValueError):reconstruct_declared_sidewall(p,f,range(8),range(6))
        with self.assertRaises(ValueError):reconstruct_declared_sidewall(p,f,range(10),[0,1,2,3,4,4])


if __name__=='__main__':unittest.main()
