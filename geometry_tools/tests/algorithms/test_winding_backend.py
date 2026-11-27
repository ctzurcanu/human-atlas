import unittest
import numpy as np
from geometry_tools.winding_backend import winding_numbers

class WindingBackendTests(unittest.TestCase):
    def source(self):
        v=np.array([[1,1,1],[-1,-1,1],[-1,1,-1],[1,-1,-1]],float)
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        return v,f
    def test_oriented_inside_and_outside_are_distinct(self):
        v,f=self.source();w=winding_numbers(v,f,np.array([[0,0,0],[4,4,4]],float))
        self.assertGreater(w[0],.99);self.assertLess(abs(w[1]),.01)
    def test_small_source_gap_does_not_destroy_the_interior(self):
        v,f=self.source();w=winding_numbers(v,f[:-1],np.array([[0,0,0],[4,4,4]],float))
        self.assertGreater(w[0],.7);self.assertLess(abs(w[1]),.1)
    def test_world_translation_is_applied_only_once(self):
        v,f=self.source();q=np.array([[0,0,0],[4,4,4]],float);offset=np.array([0,1.4,.04])
        self.assertTrue(np.allclose(winding_numbers(v,f,q),winding_numbers(v+offset,f,q+offset),atol=1e-6))

if __name__=='__main__':unittest.main()
