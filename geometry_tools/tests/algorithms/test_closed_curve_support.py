import unittest
import numpy as np
from geometry_tools.algorithms.closed_curve_support import closest_closed_curve,curve_support_planes


class ClosedCurveTests(unittest.TestCase):
    def test_annular_support_does_not_fill_the_opening_or_extrapolate_edges(self):
        curve=np.array([[-2.,-1.,0.],[2.,-1.,0.],[2.,1.,0.],[-2.,1.,0.]])
        q=np.array([[0.,0.,0.],[0.,1.,.3],[3.,2.,0.],[-2.,0.,0.]])
        r=closest_closed_curve(q,curve,batch_size=1)
        np.testing.assert_allclose(r['distance_m'],[1.,.3,np.sqrt(2.),0.])
        self.assertFalse(r['distance_m'][0] < .5)
        np.testing.assert_allclose(np.linalg.norm(q-r['points'],axis=1),r['distance_m'])
        np.testing.assert_allclose(r['points'],curve[r['edge_ids']]+r['edge_fraction'][:,None]*(np.roll(curve,-1,axis=0)-curve)[r['edge_ids']])

    def test_nonplanar_curve_uses_its_actual_edges_at_anatomical_scale(self):
        curve=np.array([[0.,0.,0.],[.004,0.,.001],[.004,.003,0.],[0.,.003,-.001]])+[0.,.92,-.03]
        midpoint=(curve[0]+curve[1])*.5
        r=closest_closed_curve(np.array([midpoint]),curve)
        self.assertLess(r['distance_m'][0],1e-16)
        self.assertAlmostEqual(r['edge_fraction'][0],.5)
        with self.assertRaises(ValueError):closest_closed_curve(np.array([midpoint]),np.r_[curve,curve[:1]])

    def test_affine_support_contains_the_nonplanar_rim_and_excludes_distant_wall(self):
        t=np.arange(20)*2*np.pi/20
        rim=np.c_[.003*np.cos(t),.002*np.sin(t),.0003*np.sin(2*t)]+[0,.92,-.03]
        center,normals,offsets=curve_support_planes(rim,.001)
        self.assertTrue(np.all(offsets-(rim-center)@normals.T>0))
        self.assertFalse(np.all(offsets-(np.array([0,.92,.03])-center)@normals.T>0))
        a,b=rim[0],rim[8];fraction=.37
        left=offsets-((1-fraction)*a+fraction*b-center)@normals.T
        right=(1-fraction)*(offsets-(a-center)@normals.T)+fraction*(offsets-(b-center)@normals.T)
        np.testing.assert_allclose(left,right,atol=2e-16)


if __name__=='__main__':unittest.main()
