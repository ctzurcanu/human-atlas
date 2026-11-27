import unittest
import numpy as np
from geometry_tools.algorithms.declared_motion_groups import declared_motion_groups
class DeclaredMotionTests(unittest.TestCase):
    def test_corresponding_gapped_profiles_keep_their_original_offsets(self):
        maps=[np.array([0,0,1]),np.array([2,2,3])]
        links=[dict(firstActor=0,firstControls=[0,1],secondActor=1,secondControls=[0,1],reason='explicit source contour')]
        result,_,proof=declared_motion_groups(maps,links,fixed_groups=[1,3])
        self.assertEqual(result[0][0],result[1][0]);self.assertNotEqual(result[0][2],result[1][2])
        points=[np.array([[0.,0,0],[0,1,0]]),np.array([[.00005,0,0],[.00005,1,0]])]
        d=np.array([.001,.002,.003]);np.testing.assert_allclose((points[1]+d)-(points[0]+d),points[1]-points[0],atol=1e-18)
        self.assertFalse(proof['sourcePointsWelded'])
    def test_a_fixed_intracranial_star_cannot_be_released_by_a_motion_link(self):
        links=[dict(firstActor=0,firstControls=[0],secondActor=1,secondControls=[0],reason='source interface')]
        with self.assertRaisesRegex(ValueError,'separately fixed'):
            declared_motion_groups([np.array([0]),np.array([1])],links,fixed_groups=[1])
if __name__=='__main__':unittest.main()
