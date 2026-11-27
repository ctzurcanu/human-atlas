import unittest
import numpy as np
from geometry_tools.algorithms.coupled_hand_cohort import coupled_hand_cohort


class CoupledHandTests(unittest.TestCase):
    def test_common_digit_translation_preserves_sections_and_fixed_palm(self):
        axis=np.array([[0.,0.,0.],[0.,0.,.03],[0.,0.,.06],[0.,0.,.08]])
        source={k:axis+[.02*k,0,0]for k in range(1,6)}
        target={k:v+[.005,0,0]for k,v in source.items()}
        ring=np.array([[.001*np.cos(t),.001*np.sin(t),.05]for t in np.linspace(0,2*np.pi,6,endpoint=False)])
        v=np.r_[np.concatenate([ring+[.02*k,0,0]for k in range(1,6)]),[[0,0,-.03]]]
        sections=[np.arange(6*k,6*k+6)for k in range(5)]
        q,proof=coupled_hand_cohort(v,sections,source,target,[0,1,0],[0,1,0],v)
        np.testing.assert_allclose(q[:-1]-v[:-1],np.tile([.005,0,0],(30,1)),atol=1e-14)
        np.testing.assert_array_equal(q[-1],v[-1])
        for s in sections:
            np.testing.assert_allclose(np.linalg.norm(q[s][:,None]-q[s][None],axis=2),np.linalg.norm(v[s][:,None]-v[s][None],axis=2),atol=1e-14)

    def test_partial_hand_cannot_be_called_a_coupled_hand(self):
        with self.assertRaises(ValueError):
            coupled_hand_cohort(np.zeros((1,3)),[],{}, {},[0,1,0],[0,1,0],np.zeros((1,3)))


if __name__=='__main__':unittest.main()
