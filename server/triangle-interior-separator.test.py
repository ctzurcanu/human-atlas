from pathlib import Path
import sys,unittest,numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ta98_source_triangle_separation import projection_separator

class TriangleInteriorSeparator(unittest.TestCase):
    def test_adjacent_coplanar_interiors_use_shared_edge_plane(self):
        a=np.array([[0,0,0],[1,0,0],[0,1,0]],float)
        b=np.array([[1,0,0],[1,1,0],[0,1,0]],float)
        n,g=projection_separator(a,b,prefer_interior_separation=True)
        self.assertAlmostEqual(g,0)
        self.assertAlmostEqual(n[2],0)
        self.assertGreater((a.mean(0)-b.mean(0))@n,.4)
        self.assertGreaterEqual((a@n).min()-(b@n).max(),-1e-14)
    def test_positive_source_gap_stays_a_valid_separator(self):
        a=np.array([[0,0,1],[1,0,1],[0,1,1]],float)
        b=a-np.array([0,0,1.])
        n,g=projection_separator(a,b,prefer_interior_separation=True)
        self.assertAlmostEqual(g,1)
        self.assertAlmostEqual((a@n).min()-(b@n).max(),g)
    def test_default_direction_is_backward_compatible(self):
        a=np.array([[0,0,0],[1,0,0],[0,1,0]],float)
        b=np.array([[1,0,0],[1,1,0],[0,1,0]],float)
        n,g=projection_separator(a,b)
        np.testing.assert_array_equal(n,[0,0,1])

if __name__=='__main__':unittest.main()
