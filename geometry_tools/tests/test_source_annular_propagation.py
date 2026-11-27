import unittest

import numpy as np

from geometry_tools.algorithms.source_annular_propagation import propagate_source_annuli
from geometry_tools.tests.test_source_profile_axis import tube
from geometry_tools.algorithms.tube_continuity import loft_faces


class SourceAnnularPropagationTests(unittest.TestCase):
    def test_recovers_full_source_strips_without_triangle_pairing(self):
        v, f = tube()
        # Reverse individual quad diagonals and permute the triangle order.
        for k in range(0, len(f), 4):
            a, c, d = f[k]
            _, _, b = f[k+1]
            f[k:k+2] = [[a, c, b], [b, c, d]]
        f = f[np.roll(np.arange(len(f)), 7)]
        before = v.copy()
        result = propagate_source_annuli(v, f, [np.arange(18, 24)])
        self.assertEqual(len(result['profiles']), 7)
        self.assertEqual(result['sourceFacesOwned'], len(f))
        faces = [i for s in result['strips'] for i in s['source_face_ids']]
        self.assertEqual(sorted(faces), list(range(len(f))))
        np.testing.assert_array_equal(v, before)

    def test_nearby_tube_stays_outside_the_declared_source_component(self):
        a, fa = tube()
        b, fb = tube(.101)
        result = propagate_source_annuli(np.r_[a, b], np.r_[fa, fb+len(a)], [np.arange(6)])
        self.assertEqual(result['sourceFacesOwned'], len(fa))
        self.assertTrue(all(np.all(r < len(a)) for r in result['profiles']))
        self.assertEqual(len(result['unassignedSourceFaceIds']), len(fb))

    def test_twelve_control_source_sections_keep_their_actual_shape(self):
        angles=np.arange(12)*2*np.pi/12
        v=np.array([[.1+.003*np.cos(t)+.0002*k*k,.001*k,.05+.001*np.sin(t)]
                    for k in range(5)for t in angles])
        f=loft_faces(5,12);f=f[np.roll(np.arange(len(f)),13)]
        before=v.copy();source_faces=f.copy()
        result=propagate_source_annuli(v,f,[np.arange(24,36)])
        self.assertEqual(len(result['profiles']),5)
        self.assertEqual(result['sourceSectionControlCounts'],[12])
        self.assertEqual(result['sourceFacesOwned'],len(f))
        self.assertTrue(all(len(r)==12 for r in result['profiles']))
        np.testing.assert_array_equal(v,before)
        np.testing.assert_array_equal(f,source_faces)

    def test_different_source_section_sizes_remain_separate_components(self):
        a,fa=tube()
        angles=np.arange(12)*2*np.pi/12
        b=np.array([[.2+.003*np.cos(t),.001*k,.05+.001*np.sin(t)]
                    for k in range(4)for t in angles])
        fb=loft_faces(4,12)+len(a)
        result=propagate_source_annuli(np.r_[a,b],np.r_[fa,fb],
                                       [np.arange(6),len(a)+np.arange(12)])
        self.assertEqual(result['sourceSectionControlCounts'],[6,12])
        self.assertEqual(result['sourceFacesOwned'],len(fa)+len(fb))
        for ring in result['profiles']:
            self.assertTrue(np.all(ring<len(a))or np.all(ring>=len(a)))


if __name__ == '__main__':
    unittest.main()
