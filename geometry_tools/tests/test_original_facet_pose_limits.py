import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.original_facet_pose_limits import OriginalFacetPoseLimits


class FacetPoseLimits(unittest.TestCase):
    def model(self):
        p=np.array([[0.,0.,0.],[.01,0.,0.],[0.,.008,.001]])
        return OriginalFacetPoseLimits(Mesh(p,np.array([[0,1,2]]),'test'),Mesh(p,np.array([[0,1,2]]),'test'),
            np.arange(3),p+np.array([.001,-.002,.001]),np.full(3,.03))

    def test_orientation_chain_matches_exact_pose_finite_differences(self):
        m=self.model();x=np.array([[.00004,-.00001,.00002],[-.00003,.00002,.003],
            [.00001,-.00003,.00004],[.0002,-.0001,.0003],[-.0001,.0003,.0002],[.0002,.0001,-.0002]])
        _,gradient=m.orientation_gradients(x,[0]);finite=np.zeros(x.size);step=1e-8
        for i in range(x.size):
            a=x.copy();b=x.copy();a.ravel()[i]+=step;b.ravel()[i]-=step
            finite[i]=(m.orientation_values(a)[0]-m.orientation_values(b)[0])/(2*step)
        np.testing.assert_allclose(gradient.toarray()[0],finite,rtol=2e-6,atol=2e-7)

    def test_actual_turn_is_rejected_and_local_row_has_correct_sign(self):
        m=self.model();x=np.zeros((6,3));x[1,2]=.003
        bad,proof=m.orientation_violations(x)
        self.assertEqual(bad.tolist(),[0]);self.assertGreater(proof['maximumFacetTurn_deg'],10.)
        A,b,_=m.local_orientation_rows(x,bad)
        self.assertGreater((b-A@x.ravel()).max(),0.)


if __name__=='__main__':unittest.main()
