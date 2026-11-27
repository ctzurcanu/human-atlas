import unittest
import numpy as np
from geometry_tools.algorithms.declared_edge_collapse import collapse_edge
from geometry_tools.algorithms.tube_continuity import continuity_audit


class DeclaredEdgeCollapseTests(unittest.TestCase):
    def source(self):
        p=np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)*1e-3
        f=np.array([[0,2,4],[2,1,4],[1,3,4],[3,0,4],
                    [2,0,5],[1,2,5],[3,1,5],[0,3,5]])
        return p,f

    def test_link_and_continuous_image_coverage(self):
        p,f=self.source();q,g,proof=collapse_edge(p,f,[0,2],owned_faces=range(len(f)),maximum_error_m=.0008)
        self.assertEqual(len(g),6)
        self.assertTrue(proof['linkConditionPassed'])
        self.assertTrue(proof['collapsedTriangleImagesCovered'])
        self.assertAlmostEqual(proof['maximumTwoWayPiecewiseLinearDistance_m'],np.sqrt(2)*.0005)
        np.testing.assert_array_equal(q[[1,3,4,5]],p[[1,3,4,5]])
        self.assertTrue(continuity_audit(q,g,expected_components=1,expected_mouths=0,expected_euler=2)['topological_continuity_pass'])

    def test_fixed_or_unowned_star_is_rejected(self):
        p,f=self.source()
        with self.assertRaises(ValueError):collapse_edge(p,f,[0,2],protected=[2],owned_faces=range(8),maximum_error_m=.0008)
        with self.assertRaises(ValueError):collapse_edge(p,f,[0,2],owned_faces=[0,1],maximum_error_m=.0008)

    def test_error_bound_cannot_expand(self):
        p,f=self.source()
        with self.assertRaises(ValueError):collapse_edge(p,f,[0,2],owned_faces=range(8),maximum_error_m=.0001)


if __name__=='__main__':unittest.main()
