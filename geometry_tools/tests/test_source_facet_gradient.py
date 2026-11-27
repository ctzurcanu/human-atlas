import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_facet_gradient import facet_gradient_rows, compact_gradient_rows, gradient_evidence


class SourceFacetGradientTests(unittest.TestCase):
    def setUp(self):
        self.mesh = Mesh(np.array([[0.,0.,0.],[2.,0.,0.],[1.,1.,0.]]), np.array([[0,1,2]]), 'test')

    def test_affine_derivative(self):
        rows = facet_gradient_rows(self.mesh, np.eye(3))
        displacement = self.mesh.vertices @ np.array([[.2,.3,.1],[-.1,.2,.4],[0,0,0]])
        np.testing.assert_allclose(rows[0] @ displacement, [[.2,.3,.1],[-.1,.2,.4]])

    def test_rigid_profile_translation_has_zero_gradient(self):
        rows = facet_gradient_rows(self.mesh, np.ones((3,1)))
        np.testing.assert_array_equal(rows, np.zeros((1,2,1)))

    def test_compaction_preserves_strongest_pair(self):
        rows = np.array([[1.,-1.,0.],[-3.,3.,0.],[0.,1.,2.],[0.,0.,0.]])
        compact = compact_gradient_rows(rows)
        self.assertEqual(len(compact), 2)
        for d in [np.eye(3), np.arange(9).reshape(3,3)]:
            self.assertAlmostEqual(np.linalg.norm(rows@d,axis=1).max(), np.linalg.norm(compact@d,axis=1).max())

    def test_small_gradient_preserves_orientation(self):
        d = np.array([[0.,0.,0.],[.2,0.,.2],[.1,-.1,.1]])
        evidence = gradient_evidence(self.mesh, np.eye(3), d)
        self.assertLess(evidence['maximumGradientFrobenius'], .5)
        old = self.mesh.triangles[0]
        new = (self.mesh.vertices+d)[self.mesh.faces[0]]
        self.assertGreater(np.cross(new[1]-new[0],new[2]-new[0]) @ np.cross(old[1]-old[0],old[2]-old[0]), 0)
