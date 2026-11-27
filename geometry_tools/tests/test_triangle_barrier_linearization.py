import unittest
import numpy as np
from scipy.sparse import csr_matrix
from geometry_tools.algorithms.triangle_barrier_linearization import oriented_triangle_barrier


class TriangleBarrierTest(unittest.TestCase):
    def test_two_profile_translation_has_exact_linear_area(self):
        v=np.array([[0.,0,0],[1,0,0],[0,1,0],[1,1,0]])
        f=np.array([[0,1,2],[1,3,2]])
        W=csr_matrix([[1.,0],[1,0],[0,1],[0,1]])
        n=np.tile([0.,0,1],(2,1));scale=np.ones(2)
        value,J=oriented_triangle_barrier(v,f,n,scale,W)
        delta=np.array([[.2,.1,.3],[-.1,.4,-.2]])
        new,_=oriented_triangle_barrier(v+W@delta,f,n,scale,W)
        np.testing.assert_allclose(new,value+J@delta.ravel(),atol=1e-14)

    def test_core_influence_derivative_matches_independent_difference(self):
        v=np.array([[.1,.2,.3],[.9,.1,.2],[.2,.8,.5]])
        W=csr_matrix([[1.,0,0],[0,1,0],[.2,.3,.5]])
        n=np.array([[.1,.2,1.]]);n/=np.linalg.norm(n)
        value,J=oriented_triangle_barrier(v,[[0,1,2]],n,[.5],W)
        delta=np.array([[.2,-.3,.1],[-.1,.2,.5],[.1,.3,-.2]])
        h=1e-7
        plus,_=oriented_triangle_barrier(v+h*(W@delta),[[0,1,2]],n,[.5],W)
        minus,_=oriented_triangle_barrier(v-h*(W@delta),[[0,1,2]],n,[.5],W)
        np.testing.assert_allclose((plus-minus)/(2*h),J@delta.ravel(),atol=1e-9)

    def test_common_translation_preserves_all_oriented_areas(self):
        v=np.array([[0.,0,0],[1,0,0],[0,1,0]])
        value,J=oriented_triangle_barrier(v,[[0,1,2]],[[0.,0,1]],[1.],csr_matrix(np.ones((3,1))))
        np.testing.assert_array_equal(J.toarray(),np.zeros((1,3)))
        self.assertEqual(value[0],1)


if __name__=='__main__':unittest.main()
