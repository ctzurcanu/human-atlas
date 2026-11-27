"""Analytical cell tests for serialization-induced boundary degeneracy."""
import unittest
import numpy as np
from geometry_tools.algorithms.encoded_boundary_chain import conform_encoded_boundary


class EncodedBoundaryTests(unittest.TestCase):
    def setUp(self):
        self.v = np.array([[0.,0.,0.], [1.,0.,0.], [0.,1.,0.], [0.,0.,1.], [.5,0.,0.]])
        self.f = np.array([[0,2,1], [0,1,3], [1,2,3], [2,0,3]])

    @staticmethod
    def volume(v, f):
        p = v[f]
        return float(np.einsum('ij,ij->i', p[:,0], np.cross(p[:,1],p[:,2])).sum()/6)

    def test_collinear_crack_restores_the_same_tetrahedron_without_a_cap(self):
        # One side subdivides the shared edge; the other does not. Serialization
        # leaves a zero-area triangle joining the two boundary subdivisions.
        raw = np.r_[self.f[[0,2,3]], [[0,4,3], [4,1,3], [0,1,4]]]
        before = self.v.copy()
        f, parents, proof = conform_encoded_boundary(self.v, raw)
        self.assertTrue(np.array_equal(self.v, before))
        self.assertEqual(len(proof['collinearBoundaryRefinements']), 1)
        self.assertEqual(proof['trueOrUnresolvedBoundaryLoops'], [])
        self.assertEqual(proof['anatomicalOpeningCapsAdded'], 0)
        self.assertEqual(self.volume(self.v, f), 1/6)
        edges = np.r_[f[:,[0,1]], f[:,[1,2]], f[:,[2,0]]]
        _, counts = np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
        self.assertTrue(np.all(counts == 2))
        self.assertTrue(np.isin(parents, np.arange(len(raw))).all())

    def test_true_open_triangle_is_preserved(self):
        raw = self.f[1:]
        f, _, proof = conform_encoded_boundary(self.v, raw)
        self.assertEqual(len(f), len(raw))
        self.assertEqual(len(proof['trueOrUnresolvedBoundaryLoops']), 1)
        self.assertEqual(proof['collinearBoundaryRefinements'], [])

    def test_opposed_internal_faces_cancel_but_same_signed_duplicates_fail(self):
        raw = np.r_[self.f, [[0,1,4], [0,4,1]], [self.f[0]], [self.f[0][[0,2,1]]]]
        f, _, proof = conform_encoded_boundary(self.v, raw)
        self.assertEqual(self.volume(self.v, f), 1/6)
        self.assertEqual(len(f), 4)
        self.assertEqual(proof['opposedCoincidentTrianglesCancelled'], 2)
        with self.assertRaisesRegex(ValueError, 'nonunit multiplicity'):
            conform_encoded_boundary(self.v, np.r_[self.f, [self.f[0]]])


if __name__ == '__main__': unittest.main()
