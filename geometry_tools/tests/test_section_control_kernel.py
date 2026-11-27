import unittest
import numpy as np
from geometry_tools.algorithms.section_control_kernel import repair_section_controls,repair_section_control_patch


class SectionControlKernelTests(unittest.TestCase):
    def test_repair_keeps_shared_controls_and_unmoved_points(self):
        p=np.array([[0.,0],[1,1],[1,0],[0,1],[1,1]])
        A=np.array([[0.],[1],[0],[0],[1]])
        q,d,r=repair_section_controls(p,[[0,1,2],[0,2,3]],A,
            minimum_double_area=.1,maximum_displacement_m=2.)
        self.assertTrue(r['passed'])
        np.testing.assert_array_equal(q,p+A@d)
        np.testing.assert_array_equal(q[1],q[4])
        np.testing.assert_array_equal(q[[0,2,3]],p[[0,2,3]])

    def test_incompatible_shared_motion_stays_failed(self):
        p=np.array([[0.,0],[1,1],[1,0],[0,1]])
        A=np.ones((4,1))
        q,d,r=repair_section_controls(p,[[0,1,2],[0,2,3]],A,
            minimum_double_area=.1,maximum_displacement_m=2.)
        self.assertFalse(r['passed'])
        np.testing.assert_array_equal(q,p)

    def test_coupled_patch_retains_shared_source_replay(self):
        p=np.array([[0.,0],[1,1],[1,0],[0,1],[1,1]])
        A=np.array([[0.],[1],[0],[0],[1]])
        q,d,r=repair_section_control_patch(p,[[0,1,2],[0,2,3]],A,
            minimum_double_area=.1,maximum_displacement_m=2.)
        self.assertTrue(r['passed'])
        np.testing.assert_array_equal(q,p+A@d)
        np.testing.assert_array_equal(q[1],q[4])

    def test_coupled_native_timeout_cannot_qualify(self):
        p=np.array([[0.,0],[1,1],[1,0]])
        with self.assertRaises(TimeoutError):
            repair_section_control_patch(p,[[0,1,2]],np.eye(3),
                minimum_double_area=.1,maximum_displacement_m=2.,wall_seconds=1e-12)


if __name__=='__main__':unittest.main()
