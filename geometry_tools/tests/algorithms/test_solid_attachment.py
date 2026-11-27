import unittest
import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.algorithms.solid_attachment import bridge_patches
from geometry_tools.audit import audit


def patch(z, reverse=False, nonplanar=False):
    p = np.array([[-.002, -.002, z], [.002, -.002, z],
                  [.002, .002, z], [-.002, .002, z]])
    if nonplanar:
        p[2, 2] += .00015
    f = np.array([[0, 1, 2], [0, 2, 3]])
    return Mesh(p, f[:, ::-1] if reverse else f, 'fixture')


def centroid(m):
    t = m.triangles
    w = np.linalg.norm(np.cross(t[:, 1]-t[:, 0], t[:, 2]-t[:, 0]), axis=1)
    return np.einsum('i,ij->j', w, t.mean(1))/w.sum()


class SolidAttachmentTests(unittest.TestCase):
    def test_finite_bridge_has_literal_opposed_contacts_and_no_self_intersection(self):
        first, second = patch(0), patch(.02, reverse=True)
        centers = np.linspace(centroid(first), centroid(second), 12)
        v, surfaces, evidence = bridge_patches(first, second, centers, np.full((12, 2), .0015))
        np.testing.assert_array_equal(v[surfaces['first_contact']], first.triangles[:, [0, 2, 1]])
        np.testing.assert_array_equal(v[surfaces['second_contact']], second.triangles[:, [0, 2, 1]])
        result = audit(Mesh(v, np.concatenate(list(surfaces.values())), 'fixture'),
                       require_closed=True, expected_components=1, expected_genus=[0])
        self.assertTrue(result['passed'], result)
        self.assertFalse(evidence['anatomy_inferred'])

    def test_nonplanar_contact_is_preserved_and_wrong_physical_side_is_rejected(self):
        first, second = patch(0, nonplanar=True), patch(.02, reverse=True)
        centers = np.linspace(centroid(first), centroid(second), 12)
        v, surfaces, _ = bridge_patches(first, second, centers, np.full((12, 2), .0015))
        np.testing.assert_array_equal(v[surfaces['first_contact']], first.triangles[:, [0, 2, 1]])
        with self.assertRaisesRegex(ValueError, 'cap faces'):
            bridge_patches(first, second, centers, np.full((12, 2), .0015), first_side=1)

    def test_turn_preserves_the_declared_elliptical_section_area(self):
        first, second = patch(0), patch(.02, reverse=True)
        t = np.linspace(0, 1, 26)
        centers = np.c_[.002*np.sin(np.pi*t)**2, np.zeros(len(t)), .02*t]
        radii = np.tile([.0012, .0008], (len(t), 1))
        v, surfaces, _ = bridge_patches(first, second, centers, radii)
        rings = v[8:].reshape(len(t)-2, 12, 3)
        np.testing.assert_allclose(rings.mean(1), centers[1:-1], atol=1e-17)
        q = rings-rings.mean(1)[:, None]
        area = np.linalg.norm(np.cross(q, np.roll(q, -1, axis=1)).sum(1), axis=1)/2
        expected = 12/2*np.sin(2*np.pi/12)*.0012*.0008
        np.testing.assert_allclose(area, expected, rtol=2e-14)
        check = audit(Mesh(v, np.concatenate(list(surfaces.values())), 'fixture'),
                      require_closed=True, expected_components=1, expected_genus=[0])
        self.assertTrue(check['passed'], check)

    def test_irregular_contact_patch_does_not_require_circular_angular_correspondence(self):
        vertices=np.array([[-.004,-.001,0],[.004,-.001,0],[.001,.002,0]])
        first=Mesh(vertices,np.array([[0,1,2]]),'fixture')
        second=Mesh(vertices+[.001,.001,.02],np.array([[0,2,1]]),'fixture')
        centers=np.linspace(centroid(first),centroid(second),16)
        v,surfaces,evidence=bridge_patches(first,second,centers,np.tile([.002,.001],(16,1)),section_controls=16,projected_phase_alignment=False)
        np.testing.assert_array_equal(v[surfaces['first_contact']],first.triangles[:,[0,2,1]])
        np.testing.assert_array_equal(v[surfaces['second_contact']],second.triangles[:,[0,2,1]])
        result=audit(Mesh(v,np.concatenate(list(surfaces.values())),'fixture'),require_closed=True,expected_components=1,expected_genus=[0])
        self.assertTrue(result['passed'],result)
        self.assertFalse(evidence['last_strip']['projected_phase_alignment'])


if __name__ == '__main__':
    unittest.main()
