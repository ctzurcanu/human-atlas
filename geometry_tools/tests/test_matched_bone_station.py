import unittest
import numpy as np
from geometry_tools.commands.pose_hand_vessel_courses import trim_to_matched_bone_station
from geometry_tools.algorithms.source_profile_axis import complete_source_sections


class MatchedStationTests(unittest.TestCase):
    def test_unpaired_source_port_is_retained_as_one_rigid_ring(self):
        g=dict(profiles=[np.arange(6)],boundary_raw=[np.arange(6),np.arange(6,12)])
        result=complete_source_sections(g,np.arange(12))
        self.assertEqual(len(result),2)
        np.testing.assert_array_equal(result[1],np.arange(6,12))

    def test_partly_overlapping_rings_move_as_one_shared_closure(self):
        g=dict(profiles=[np.arange(6)],boundary_raw=[np.arange(3,9)])
        result=complete_source_sections(g,np.arange(9))
        self.assertEqual(len(result),1)
        np.testing.assert_array_equal(result[0],np.arange(9))

    def test_long_reference_prefix_is_trimmed_to_actual_source_station(self):
        source=np.array([[0.,0.,0.],[0.,.02,0.],[0.,.04,0.]])
        target=source*1.2+np.array([.1,.2,.3])
        guide=np.array([[.103,.19,.3],[.103,.212,.3],[.103,.236,.3],[.103,.26,.3]])
        trimmed,proof=trim_to_matched_bone_station(guide,np.array([.002,.03,0.]),source,target)
        np.testing.assert_allclose(trimmed[0],[.103,.236,.3],atol=1e-14)
        np.testing.assert_array_equal(trimmed[-1],guide[-1])
        self.assertEqual(proof['sourceBoneSegment'],1)
        self.assertAlmostEqual(proof['segmentFraction'],.5)

    def test_missing_station_does_not_guess_a_mouth(self):
        bones=np.array([[0.,0.,0.],[0.,.02,0.]])
        with self.assertRaisesRegex(ValueError,'misses'):
            trim_to_matched_bone_station(np.array([[.001,.03,0.],[.001,.04,0.]]),
                                        np.array([0.,.01,0.]),bones,bones)


if __name__=='__main__':unittest.main()
