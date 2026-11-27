import unittest
import numpy as np
from geometry_tools.algorithms.interpolating_rigid_field import interpolate_rigid_field


class InterpolatingRigidFieldTests(unittest.TestCase):
    def controls(self):
        return np.array([[0,0,0],[.001,0,0],[0,.001,0],[0,0,.001],[.001,.001,.001]],float)

    def test_interpolation_is_literal_at_all_declared_controls(self):
        p=self.controls();q=p.copy();q[1]+=[0,.0001,0]
        result,proof=interpolate_rigid_field(p,p,q,neighbours=4)
        np.testing.assert_array_equal(result,q)
        self.assertEqual(proof['exactInterpolatedControls'],len(p))

    def test_a_rigid_pose_is_recovered_away_from_controls(self):
        p=self.controls();R=np.array([[0,-1,0],[1,0,0],[0,0,1.]])
        q=p@R.T+[.3,.8,.04];x=np.array([[.0002,.0003,.0004],[.002,.001,.003]])
        result,_=interpolate_rigid_field(x,p,q,neighbours=5,batch_size=1)
        np.testing.assert_allclose(result,x@R.T+[.3,.8,.04],atol=1e-14)

    def test_conflicting_source_aliases_are_rejected(self):
        p=self.controls();q=np.r_[p,p[:1]];target=q.copy();target[-1,0]+=.001
        with self.assertRaisesRegex(ValueError,'contradictory'):
            interpolate_rigid_field(p,q,target,neighbours=4)

    def test_explicit_distant_support_remains_fixed(self):
        p=self.controls();x=np.array([[1.,1.,1.]])
        result,_=interpolate_rigid_field(x,p,p+.01,neighbours=4,support_m=(.005,.02))
        np.testing.assert_array_equal(result,x)


if __name__=='__main__':unittest.main()
