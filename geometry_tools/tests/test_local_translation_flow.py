import unittest
import numpy as np
from geometry_tools.algorithms.local_translation_flow import carry_points

class TranslationFlowTests(unittest.TestCase):
    def test_station_preserves_remote_attachment_and_variational_derivative(self):
        p=np.array([[.01,.91,.05],[.02,.928,.05],[.02,.94,.05]])
        centres=np.array([[.01,.91,.05]])
        deltas=np.array([[.005,.006,.004]])
        options=dict(upper_station=(1,.921,.934),wall_seconds=2.)
        q,J=carry_points(p,centres,deltas,**options)
        np.testing.assert_array_equal(q[-1],p[-1])
        np.testing.assert_array_equal(J[-1],np.eye(3))
        for axis in range(3):
            step=np.eye(3)[axis]*1e-7
            plus,_=carry_points(p+step,centres,deltas,**options)
            minus,_=carry_points(p-step,centres,deltas,**options)
            np.testing.assert_allclose((plus-minus)/2e-7,J[:,:,axis],rtol=1e-5,atol=1e-7)
        self.assertTrue(np.all(np.linalg.det(J)>0))

    def test_core_translation_retains_calibre(self):
        p=np.array([[.001,0,0],[-.001,0,0],[0,.001,0]])
        d=np.array([[.004,.003,.002]])
        q,J=carry_points(p,np.zeros((1,3)),d)
        np.testing.assert_allclose(q,p+d,atol=1e-15)
        np.testing.assert_allclose(J,np.broadcast_to(np.eye(3),J.shape),atol=1e-15)

if __name__=='__main__':unittest.main()
