"""Independent continuous-depth checks of the optional repair constraint mode."""
from pathlib import Path
import sys, unittest, numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ta98_station_prism_constraints import StationPrismConstraints
from ta98_shared_wall_layers import prism_jacobian_maxima

class PrismRepair(unittest.TestCase):
    def setUp(self):
        self.outer=np.array([[1,0,0],[1,1,0],[1,0,1],[1,1,1]],float)
        self.faces=np.array([[0,1,2],[1,3,2]])
        self.station=np.array([0,1,0,1]);self.fixed=np.array([True,False])

    def test_default_still_rejects_invalid_source(self):
        inner=self.outer.copy();inner[:,0]+=.1
        with self.assertRaises(AssertionError):
            StationPrismConstraints(np.array([self.outer,inner]),self.faces,self.station,self.fixed)

    def test_repair_targets_are_strictly_negative(self):
        inner=self.outer.copy();inner[:,0]-=.1;inner[1,0]=1.03
        q=StationPrismConstraints(np.array([self.outer,inner]),self.faces,self.station,self.fixed,repair_invalid=True)
        self.assertTrue(np.all(q.limit<0))
        self.assertTrue(np.any(q.evaluate(np.zeros((2,3)))[0]>0))

    def test_analytic_maxima_match_explicit_translated_geometry(self):
        inner=self.outer.copy();inner[:,0]-=.1;inner[1,0]=1.03
        fields=np.array([self.outer,inner])
        q=StationPrismConstraints(fields,self.faces,self.station,self.fixed,repair_invalid=True)
        translation=np.array([[0,0,0],[.012,-.031,.009]])
        values,_=q.evaluate(translation)
        moved=fields+translation[self.station][None]
        np.testing.assert_allclose(values[0],prism_jacobian_maxima(moved[0],moved[1],self.faces),atol=1e-14,rtol=1e-14)

if __name__=='__main__':unittest.main()
