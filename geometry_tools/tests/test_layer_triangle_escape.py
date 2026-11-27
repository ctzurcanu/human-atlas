import unittest
import numpy as np
from geometry_tools.algorithms.layer_triangle_escape import layer_triangle_exit_axes


class LayerTriangleEscape(unittest.TestCase):
    def test_oblique_faces_do_not_require_full_axis_projection_width(self):
        a=np.array([[[-.001,-.01,0.],[.001,.01,0.],[.002,.01,0.]]])
        b=np.array([[[-.002,-.01,0.],[-.001,-.01,0.],[.001,.01,0.]]])+[.00005,0.,0.]
        direction=np.array([[1.,0.,0.]])
        axis,amount,proof=layer_triangle_exit_axes(a,b,[[0.,0.,1.]],direction,clearance_m=.00002)
        old=.00002-(a[:,:,0].min()-b[:,:,0].max())
        self.assertLess(amount[0],old*.1)
        moved=a+direction[:,None]*amount[:,None,None]
        gap=np.einsum('pvi,pi->pv',moved,axis).min(1)-np.einsum('pvi,pi->pv',b,axis).max(1)
        self.assertGreaterEqual(gap[0],.00002-1e-12)
        self.assertLess(proof['maximumLayerNormalComponent'],1e-14)

    def test_normal_only_escape_is_rejected(self):
        t=np.array([[[0.,0.,0.],[.01,0.,0.],[0.,.01,0.]]])
        with self.assertRaises(ValueError):
            layer_triangle_exit_axes(t,t,[[0.,0.,1.]],[[0.,0.,1.]])


if __name__=='__main__':unittest.main()
