import unittest,time
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.model import Mesh
from geometry_tools.algorithms.fixed_interval_rows import preservation_rows
from geometry_tools.algorithms.regional_pair_intervals import regional_pair_preservation
from geometry_tools.algorithms.section_pair_separation import separated_triangle_axes


class RegionalPairIntervals(unittest.TestCase):
    def pair(self):
        p=np.array([[0.,0.,0.],[0.,.01,0.],[0.,0.,.01]])
        return Mesh(p+[.0001,0,0],np.array([[0,1,2]]),'test'),Mesh(p,np.array([[0,1,2]]),'test')

    def test_common_rotation_can_break_old_plane_while_actual_gap_is_preserved(self):
        a,b=self.pair();g=np.arange(3);h=np.arange(3,6)
        _,_,ledger,proof=regional_pair_preservation(a,b,g,h,6,
            first_region=np.ones(3,bool),second_region=np.ones(3,bool),maximum_motion_m=.001)
        self.assertEqual(proof['replacedLocalPositivePairSurrogates'],1)
        self.assertTrue(proof['actualEveryOriginalPairVerificationRequired'])
        self.assertFalse(proof['clearanceEstablished'])
        R=Rotation.from_euler('z',2,degrees=True).as_matrix()
        p=a.triangles@R.T;q=b.triangles@R.T
        self.assertLess(p[:,:,0].min()-q[:,:,0].max(),0.)
        _,gap=separated_triangle_axes(p,q)
        self.assertAlmostEqual(gap[0],.0001,places=14)
        self.assertEqual(ledger['pairs'].tolist(),[[0,0]])

    def test_nonlocal_and_fixed_controls_keep_literal_original_constraints(self):
        a,b=self.pair();g=np.arange(3);h=np.arange(3,6)
        for first_region,groups in [(np.zeros(3,bool),g),(np.ones(3,bool),np.array([-1,1,2]))]:
            old,b0,_=preservation_rows(a,b,groups,h,6,maximum_motion_m=.001)
            new,b1,ledger,_=regional_pair_preservation(a,b,groups,h,6,
                first_region=first_region,second_region=np.ones(3,bool),maximum_motion_m=.001)
            self.assertEqual((old-new).nnz,0);np.testing.assert_array_equal(b0,b1)
            self.assertEqual(len(ledger['pairs']),0)


if __name__=='__main__':unittest.main()
