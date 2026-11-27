import unittest
import numpy as np
from geometry_tools.algorithms.recipient_clearance_field import RecipientClearanceField
from geometry_tools.tests.test_sparse_group_clearance import cube


class RecipientClearanceFieldTests(unittest.TestCase):
    def test_overlap_of_two_closed_recipients_is_still_inside_their_union(self):
        field=RecipientClearanceField([(cube([0,0,0],[.01,.01,.01]),True,'a'),
                                      (cube([.005,0,0],[.015,.01,.01]),True,'b')])
        q=np.array([[.002,.005,.005],[.0075,.005,.005],[.013,.005,.005],[.02,.005,.005]])
        gap=field.clearance(q)
        self.assertTrue(np.all(gap[:3]<0));self.assertGreater(gap[-1],.0049)
        self.assertGreaterEqual(field.occupancy_points,4)

    def test_segments_are_checked_between_clear_endpoints(self):
        field=RecipientClearanceField([(cube([0,0,0],[.01,.01,.01]),True,'fixed')])
        self.assertFalse(field.segment_clear(np.array([[-.01,.005,.005]]),np.array([[.02,.005,.005]]))[0])
        self.assertTrue(field.segment_clear(np.array([[-.01,.02,.005]]),np.array([[.02,.02,.005]]))[0])


if __name__=='__main__':unittest.main()
