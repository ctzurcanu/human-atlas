import unittest
import numpy as np
from geometry_tools.algorithms.quantization_swept_volume import swept_volume_bound

class SweptVolume(unittest.TestCase):
    def test_triangle_sweep_contains_known_prism(self):
        v=np.array([[0.,0.,0.],[2.,0.,0.],[0.,3.,0.]])
        r=swept_volume_bound(v,v+[0.,0.,.1],np.array([[0,1,2]]))
        self.assertGreaterEqual(r['sweptVolumeUpperBoundMetersCubed'],.3)
        self.assertAlmostEqual(r['sweptVolumeUpperBoundMetersCubed'],.6,places=10)
    def test_oblique_sweep_and_translation_invariance(self):
        v=np.array([[0.,0.,0.],[.002,.001,0.],[0.,.001,.003]])
        b=v+np.array([[0.,2e-9,0.],[3e-9,0.,0.],[0.,0.,-2e-9]])
        f=np.array([[0,1,2]])
        a=swept_volume_bound(v,b,f);r=swept_volume_bound(v+1.25,b+1.25,f)
        self.assertGreater(a['sweptVolumeUpperBoundMetersCubed'],0)
        self.assertAlmostEqual(a['maximumControlMotionMeters'],3e-9,places=18)
        self.assertLess(abs(a['sweptVolumeUpperBoundMetersCubed']-r['sweptVolumeUpperBoundMetersCubed'])/a['sweptVolumeUpperBoundMetersCubed'],1e-6)
    def test_invalid_input(self):
        with self.assertRaises(ValueError):swept_volume_bound(np.zeros((3,3)),np.zeros((3,3)),np.array([[0,1,2]]))
if __name__=='__main__':unittest.main()
