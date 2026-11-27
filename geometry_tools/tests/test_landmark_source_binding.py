import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.landmark_source_binding import verify_landmark_source


class LandmarkSourceBindingTests(unittest.TestCase):
    def mesh(self, shift=0., frame='test:m'):
        return Mesh(np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])+shift,
                    np.array([[0,1,2]]), frame)

    def test_exact_original_source_passes(self):
        self.assertTrue(verify_landmark_source({'artery':self.mesh()}, {'artery':self.mesh()})['complete'])

    def test_same_frame_with_previous_pose_is_rejected(self):
        with self.assertRaisesRegex(ValueError, 'pose mismatch'):
            verify_landmark_source({'artery':self.mesh(.02)}, {'artery':self.mesh()})

    def test_export_tolerance_must_be_explicit(self):
        with self.assertRaises(ValueError):
            verify_landmark_source({'artery':self.mesh(1e-7)}, {'artery':self.mesh()})
        verify_landmark_source({'artery':self.mesh(1e-7)}, {'artery':self.mesh()}, coordinate_tolerance_m=1e-6)
        with self.assertRaises(ValueError):
            verify_landmark_source({}, {}, coordinate_tolerance_m=.001)

    def test_changed_ancestry_or_frame_is_rejected(self):
        changed=self.mesh()
        changed.faces=np.array([[0,2,1]])
        with self.assertRaisesRegex(ValueError,'ancestry'):
            verify_landmark_source({'artery':changed},{'artery':self.mesh()})
        with self.assertRaisesRegex(ValueError,'frame'):
            verify_landmark_source({'artery':self.mesh(frame='other:m')},{'artery':self.mesh()})


if __name__=='__main__':
    unittest.main()
