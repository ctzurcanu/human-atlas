import unittest
import numpy as np
from geometry_tools.algorithms.geodesic_contour_axis import geodesic_contour_axis


def tube(centres,radius=.0006):
    centres=np.asarray(centres,float);angle=np.arange(12)*2*np.pi/12
    v=(centres[:,None]+np.c_[radius*np.cos(angle),np.zeros(12),radius*np.sin(angle)][None]).reshape(-1,3)
    f=[]
    for i in range(len(centres)-1):
        for j in range(12):
            a=i*12+j;b=i*12+(j+1)%12;c=a+12;d=b+12
            f.extend([[a,b,d],[a,d,c]] if (i+j)%2 else [[a,b,c],[b,d,c]])
    return v,np.asarray(f,int)


class GeodesicContourAxisTests(unittest.TestCase):
    def test_curved_trunk_follows_its_actual_surface_with_irregular_diagonals(self):
        y=np.linspace(0,.04,41);centres=np.c_[.004*np.sin(np.pi*y/.04),y,np.zeros(len(y))]
        v,f=tube(centres);before=v.copy();faces=f.copy()
        axis,proof=geodesic_contour_axis(v,f,np.arange(12),centres[-1],spacing_m=.001)
        # The three-point chord misses this curved trunk by millimetres. Its
        # contour axis must follow the real source, independent of diagonals.
        expected=.004*np.sin(np.pi*axis[:,1]/.04)
        self.assertLess(np.max(abs(axis[:,0]-expected)),.00015)
        self.assertLess(np.max(abs(axis[:,2])),.00015)
        np.testing.assert_array_equal(v,before);np.testing.assert_array_equal(f,faces)
        self.assertEqual(proof['sourceBranchesRemoved'],0)
        self.assertTrue(all(len(r['sourceFaces'])>=12 for r in proof['contourWitnesses']))

    def test_a_disconnected_other_tube_is_rejected_instead_of_joined(self):
        v,f=tube([[0,0,0],[0,.01,0],[0,.02,0]])
        with self.assertRaisesRegex(ValueError,'connected physical component'):
            geodesic_contour_axis(np.r_[v,v+[.01,0,0]],np.r_[f,f+len(v)],np.arange(12),[0,.02,0])


if __name__=='__main__':unittest.main()
