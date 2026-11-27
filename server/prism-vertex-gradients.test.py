from pathlib import Path
import sys,unittest,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ta98_prism_vertex_gradients import prism_vertex_gradients
from ta98_shared_wall_layers import prism_jacobian_maxima

class PrismGradients(unittest.TestCase):
    def setUp(self):
        rng=np.random.default_rng(638)
        self.fields=rng.normal(size=(4,7,3))
        self.faces=np.array([[0,1,2],[2,3,4],[4,5,6]])
    def test_matches_independent_extrema(self):
        value,_,_=prism_vertex_gradients(self.fields,self.faces)
        expected=np.array([prism_jacobian_maxima(a,b,self.faces) for a,b in zip(self.fields[:-1],self.fields[1:])])
        np.testing.assert_allclose(value,expected,atol=1e-12,rtol=1e-12)
    def test_vertex_derivatives_match_central_differences(self):
        value,ga,gb=prism_vertex_gradients(self.fields,self.faces)
        rng=np.random.default_rng(92);direction=rng.normal(size=self.fields.shape);epsilon=1e-6
        numerical=(prism_vertex_gradients(self.fields+epsilon*direction,self.faces)[0]-prism_vertex_gradients(self.fields-epsilon*direction,self.faces)[0])/(2*epsilon)
        exact=np.einsum('kfvi,kfvi->kf',ga,direction[:-1,self.faces])+np.einsum('kfvi,kfvi->kf',gb,direction[1:,self.faces])
        np.testing.assert_allclose(exact,numerical,atol=1e-7,rtol=1e-7)
    def test_common_translation_has_zero_derivative(self):
        _,ga,gb=prism_vertex_gradients(self.fields,self.faces)
        np.testing.assert_allclose((ga+gb).sum(2),0,atol=1e-14)

if __name__=='__main__':unittest.main()
