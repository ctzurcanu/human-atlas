import unittest
import numpy as np
from geometry_tools.algorithms.skin_transition_guard import skin_transition_guard


class SkinTransitionGuardTests(unittest.TestCase):
    def gate(self,old,new,op,npairs,known=None):
        k=np.ones(len(old),bool)
        return skin_transition_guard(old,k,new,k if known is None else known,op,npairs)

    def test_partial_improvement_retains_unresolved_exit(self):
        p=self.gate([.001,.002,-.001],[-.001,.002,-.001],[[0,1],[1,2]],[[1,2]])
        self.assertTrue(p['passed']);self.assertEqual(p['afterExteriorControls'],1)

    def test_new_triangle_crossing_rejects_fewer_exterior_controls(self):
        self.assertFalse(self.gate([.001,.002],[-.001,-.001],[[1,2]],[[2,3]])['passed'])

    def test_existing_exit_cannot_get_worse(self):
        self.assertFalse(self.gate([.001,.002],[-.001,.003],[[1,2]],[])['passed'])

    def test_unknown_classification_cannot_hide_a_defect(self):
        self.assertFalse(self.gate([.001,.002],[-.001,-.001],[[1,2]],[],np.array([True,False]))['passed'])

    def test_identity_is_not_an_improvement(self):
        self.assertFalse(self.gate([.001],[.001],[[1,2]],[[1,2]])['passed'])
