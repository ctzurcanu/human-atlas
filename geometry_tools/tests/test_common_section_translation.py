import unittest
import numpy as np
import shapely
from geometry_tools.algorithms.common_section_translation import common_section_translations


class CommonTranslationTests(unittest.TestCase):
    def test_same_offset_must_clear_every_member(self):
        square=np.array([[0.,0],[1,0],[1,1],[0,1]])
        context=dict(container=shapely.box(-1,-1,8,3),solid=shapely.box(2.5,-.5,3.5,1.5))
        valid,r=common_section_translations([square,square+[2,0]],[[0,0],[4,0]],context,
            clearance_m=.01,maximum_motion_m=5)
        np.testing.assert_array_equal(valid,[1])
        self.assertTrue(r['relativeSourceGeometryLiteral'])

    def test_source_crossing_is_not_fixed_by_common_translation(self):
        square=np.array([[0.,0],[1,0],[1,1],[0,1]])
        valid,r=common_section_translations([square,square+[.5,0]],[[0,0]],
            dict(container=shapely.box(-2,-2,4,4)),clearance_m=.01)
        self.assertEqual(len(valid),0)
        self.assertFalse(r['sourceMutualClearance'])

    def test_declared_source_gap_preservation_does_not_enlarge_close_peers(self):
        square=np.array([[0.,0],[1,0],[1,1],[0,1]])
        valid,r=common_section_translations([square,square+[1.005,0]],[[0,0]],
            dict(container=shapely.box(-2,-2,4,4)),clearance_m=.01,
            preserve_source_mutual_gaps=True)
        np.testing.assert_array_equal(valid,[0])
        self.assertAlmostEqual(r['sourceMutualGaps_m'][0],.005)


if __name__=='__main__':unittest.main()
