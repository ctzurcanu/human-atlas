import unittest
import numpy as np
from geometry_tools.algorithms.contact_profile_cohort import contact_profile_cohort
from geometry_tools.distance import triangle_pairs


class ContactCohortTests(unittest.TestCase):
    def test_intersection_survives_common_profile_transport(self):
        a=np.array([[0.,0,0],[1,0,0],[0,1,0]])
        b=np.array([[.25,.25,-1],[.25,.25,1],[.75,.25,0]])
        r=contact_profile_cohort([3,3],[np.array([[0,1,2]])]*2,[(0,1,[[0,0]])])
        self.assertEqual(r['clusters'],1)
        delta=np.array([.2,-.4,.1]); before=triangle_pairs(a[None],b[None])
        self.assertEqual(before['relation'].tolist(),['crossing'])
        after=triangle_pairs((a+delta)[None],(b+delta)[None])
        self.assertEqual(before['relation'].tolist(),after['relation'].tolist())
        self.assertFalse(r['sourceMeshesMerged'])

    def test_external_fixed_profile_protects_the_complete_contact(self):
        r=contact_profile_cohort([4,4],[np.array([[0,1,2]])]*2,
                                 [(0,1,[[0,0]])],fixed=[(0,1)])
        c=r['clusterOfGroup']; self.assertEqual(c[0][:3],c[1][:3])
        self.assertIn(c[1][2],r['fixedClusters'])
        self.assertNotIn(c[1][3],r['fixedClusters'])

    def test_exact_aliases_couple_without_contact_or_tolerance_welding(self):
        r=contact_profile_cohort([2,2],[np.empty((0,3),int)]*2,[],
                                aliases=[[(0,1),(1,0)]])
        c=r['clusterOfGroup'];self.assertEqual(c[0][1],c[1][0])
        self.assertNotEqual(c[0][0],c[1][1])


if __name__=='__main__':unittest.main()
