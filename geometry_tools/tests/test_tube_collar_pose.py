import unittest
import numpy as np
from geometry_tools.algorithms.tube_collar_pose import pose_parent_collar


class TubeCollarTests(unittest.TestCase):
    def ring(self,z=0.):
        t=np.arange(6)*np.pi/3
        return np.c_[.001*np.cos(t),.001*np.sin(t),np.full(6,z)]

    def test_collar_is_opposite_child_and_retains_calibre(self):
        mouth=self.ring();source=self.ring(-.002)
        target=mouth+np.array([.3,.8,.04]);child=target+np.array([.001,0.,.002])
        result,proof=pose_parent_collar(mouth,source,target,child)
        self.assertAlmostEqual(proof['parentChildDirectionDot'],-1.)
        self.assertAlmostEqual(np.linalg.norm(result.mean(0)-target.mean(0)),.002)
        np.testing.assert_allclose(result-result.mean(0),source-source.mean(0),atol=1e-15)

    def test_distorted_mouth_cannot_be_hidden_by_collar(self):
        a=self.ring();b=a.copy();b[0,0]+=.00005
        with self.assertRaisesRegex(ValueError,'mouth shape'):
            pose_parent_collar(a,self.ring(-.002),b,self.ring(.002))

    def test_missing_child_interval_is_rejected(self):
        a=self.ring()
        with self.assertRaisesRegex(ValueError,'nonzero'):
            pose_parent_collar(a,self.ring(-.002),a,a)


if __name__=='__main__':unittest.main()
