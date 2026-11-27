import unittest
import numpy as np
from geometry_tools.algorithms.polyline_distance import directed_polyline_distance, polyline_hausdorff


class CompletePolylineDistanceTests(unittest.TestCase):
    def test_same_contour_different_sampling(self):
        a = [[0,0,0],[2,0,0],[2,2,0],[0,2,0]]
        b = [[0,0,0],[.2,0,0],[2,0,0],[2,1,0],[2,2,0],[0,2,0]]
        r = polyline_hausdorff(a,b,tolerance=1e-9)
        self.assertLess(r["upper"],1e-9)

    def test_maximum_can_be_inside_source_edge(self):
        # Endpoints lie on the target V; the middle of the source is farther.
        r = directed_polyline_distance([[-1,0,0],[1,0,0]],
            [[-1,0,0],[0,1,0],[1,0,0]], source_closed=False,
            target_closed=False,tolerance=1e-8)
        exact = 1/np.sqrt(2)
        self.assertLessEqual(r["lower"],exact+1e-14)
        self.assertGreaterEqual(r["upper"],exact-1e-14)
        self.assertLessEqual(r["upper"]-r["lower"],1e-8)
        np.testing.assert_allclose(r["witness"]["sourcePoint"],[0,0,0])

    def test_degenerate_edges_and_parallel_distance(self):
        r = polyline_hausdorff([[0,0,0],[0,0,0],[1,0,0]],
                              [[0,.3,0],[1,.3,0]],
                              source_closed=False,target_closed=False,tolerance=1e-10)
        self.assertAlmostEqual(r["lower"],.3)
        self.assertLess(r["upper"]-.3,1e-10)

    def test_limit_failure_produces_no_uncertified_result(self):
        with self.assertRaises(RuntimeError):
            directed_polyline_distance([[-1,0,0],[1,0,0]],
                [[-1,0,0],[0,1,0],[1,0,0]], source_closed=False,
                target_closed=False,max_evaluations=2)

    def test_nonfinite_tolerance_cannot_claim_a_certificate(self):
        with self.assertRaises(ValueError):
            directed_polyline_distance([[0,0],[1,0]], [[0,1],[1,1]], tolerance=float("nan"))


if __name__ == "__main__":
    unittest.main()
