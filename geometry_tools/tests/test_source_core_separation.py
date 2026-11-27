import unittest
import numpy as np
from geometry_tools.algorithms.source_core_separation import core_separation_candidates
from geometry_tools.distance import triangle_pairs


class CoreSeparationTests(unittest.TestCase):
    def fixture(self):
        v=np.array([[-1,-1,0],[1,-1,0],[0,1,0],
                    [0,-.5,-1],[0,-.5,1],[0,.5,0],[8,8,8]],float)*.001
        return v,np.array([[0,1,2],[3,4,5]])

    def test_source_triangle_plane_has_real_gap_and_fixed_recipient(self):
        v,f=self.fixture()
        self.assertEqual(triangle_pairs(v[f[:1]],v[f[1:]])['relation'][0],'crossing')
        candidates=core_separation_candidates(v,f,[0,1],[0,1,2],gap_m=.000005,maximum_motion_m=.002)
        self.assertTrue(candidates)
        for q,p in candidates:
            np.testing.assert_array_equal(q[3:],v[3:])
            n=np.array(p['axis']);gap=(q[f[0]]@n).min()-(q[f[1]]@n).max()
            self.assertGreaterEqual(gap,.000005-1e-10)
            self.assertNotIn(triangle_pairs(q[f[:1]],q[f[1:]])['relation'][0],['crossing','coplanar-overlap'])
            self.assertLessEqual(np.linalg.norm(q-v,axis=1).max(),.002)

    def test_unrelated_free_control_cannot_move_fixed_crossing_faces(self):
        v,f=self.fixture()
        self.assertEqual(core_separation_candidates(v,f,[0,1],[6],gap_m=.000005,maximum_motion_m=.002),[])

    def test_fixed_recipient_barrier_preserves_the_original_interval(self):
        v,f=self.fixture();target=v[f[0]]-[0,0,.00001]
        barrier=dict(face=0,triangle=target,axis=[0,0,1],minimum_gap_m=.00001)
        candidates=core_separation_candidates(v,f,[0,1],[0,1,2],gap_m=.000005,maximum_motion_m=.002,barriers=[barrier])
        self.assertTrue(candidates)
        for q,_ in candidates:
            np.testing.assert_array_equal(q[3:],v[3:])
            self.assertGreaterEqual(float(q[:3,2].min()-target[:,2].max()),.00001-1e-10)

    def test_an_already_crossed_or_unnormalized_barrier_is_rejected(self):
        v,f=self.fixture()
        for axis,gap in [([0,0,1],.00002),([0,0,2],.00001)]:
            with self.subTest(axis=axis,gap=gap),self.assertRaises(ValueError):
                core_separation_candidates(v,f,[0,1],[0,1,2],gap_m=.000005,maximum_motion_m=.002,
                    barriers=[dict(face=0,triangle=v[f[0]]-[0,0,.00001],axis=axis,minimum_gap_m=gap)])

    def test_rejects_invalid_indices_and_unbounded_or_nonfinite_limits(self):
        v,f=self.fixture()
        for pair,movable,gap,maximum,seconds in [([0,2],[0],1e-6,1e-4,3),([0,1],[-1],1e-6,1e-4,3),
                ([0,1],[0,0],1e-6,1e-4,3),([0,1],[0],1e-6,float('inf'),3),
                ([0,1],[0],1e-6,1e-4,6),([0,1],[0],-1e-6,1e-4,3)]:
            with self.subTest(pair=pair,movable=movable,maximum=maximum,seconds=seconds),self.assertRaises(ValueError):
                core_separation_candidates(v,f,pair,movable,gap_m=gap,maximum_motion_m=maximum,wall_seconds=seconds)
