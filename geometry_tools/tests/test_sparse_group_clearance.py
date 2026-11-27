import unittest

import numpy as np

from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.sparse_group_clearance import sparse_group_clearance


def cube(low, high):
    a, b, c = low; x, y, z = high
    v = np.array([[a,b,c], [x,b,c], [x,y,c], [a,y,c],
                  [a,b,z], [x,b,z], [x,y,z], [a,y,z]])
    f = np.array([[0,2,1], [0,3,2], [4,5,6], [4,6,7],
                  [0,1,5], [0,5,4], [1,2,6], [1,6,5],
                  [2,3,7], [2,7,6], [3,0,4], [3,4,7]])
    return Mesh(v, f, 'test:metres')


class SparseClearanceTests(unittest.TestCase):
    def test_clearance_preserves_rigid_polygon(self):
        bone = SurfaceIndex(cube([0.,0.,0.], [.01,.01,.01]))
        skin = SurfaceIndex(cube([-.02,-.02,-.02], [.03,.03,.03]))
        v = np.array([[.0101,.004,.004], [.0101,.006,.004], [.0101,.005,.006]])
        q, proof = sparse_group_clearance(v, np.array([[0,1,2]]), [np.arange(3)], bone, skin)
        self.assertTrue(proof['complete'], proof)
        self.assertGreaterEqual(proof['minimumSampledBoneGap_m'], .0002-1e-8)
        np.testing.assert_allclose(q-v, np.broadcast_to(q[0]-v[0], q.shape), atol=1e-15)
        np.testing.assert_allclose(np.linalg.norm(q[:,None]-q[None],axis=2),
                                   np.linalg.norm(v[:,None]-v[None],axis=2), atol=1e-15)

    def test_impossible_fixed_clearance_retains_source(self):
        shift = np.array([0.,1.1,0.])
        b = cube([0.,0.,0.], [.01,.01,.01]); bone = Mesh(b.vertices+shift, b.faces, b.frame)
        s = cube([-.02,-.02,-.02], [.03,.03,.03]); skin = Mesh(s.vertices+shift, s.faces, s.frame)
        v = np.array([[.0101,.004,.004], [.0101,.006,.004], [.0101,.005,.006]])+shift
        q, proof = sparse_group_clearance(v, np.array([[0,1,2]]), [np.arange(3)],
                                           SurfaceIndex(bone), SurfaceIndex(skin))
        self.assertFalse(proof['complete'])
        np.testing.assert_array_equal(q, v)

    def test_crossed_rigid_section_can_take_one_exterior_direction(self):
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        skin=SurfaceIndex(cube([-.02,-.02,-.02],[.03,.03,.03]))
        v=np.array([[.009,.004,.004],[.011,.006,.004],[.011,.005,.006]])
        q,proof=sparse_group_clearance(v,np.array([[0,1,2]]),[np.arange(3)],bone,skin,
                                       coherent_bone_sections=True)
        self.assertTrue(proof['complete'],proof)
        self.assertGreaterEqual(q[:,0].min(),.0102-1e-8)
        np.testing.assert_allclose(q-v,np.broadcast_to(q[0]-v[0],q.shape),atol=1e-15)

    def test_unrelated_fixed_component_is_not_moved_or_claimed_clear(self):
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        skin=SurfaceIndex(cube([-.02,-.02,-.02],[.03,.03,.03]))
        moving=np.array([[.0101,.004,.004],[.0101,.006,.004],[.0101,.005,.006]])
        fixed=moving-[.005,0,0];v=np.r_[moving,fixed]
        q,proof=sparse_group_clearance(v,np.array([[0,1,2],[3,4,5]]),[np.arange(3),np.arange(3,6)],bone,skin,active_raw=np.arange(3))
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_array_equal(q[3:],v[3:])
        self.assertGreater(proof['uncheckedFixedSamples'],0)


if __name__ == '__main__':
    unittest.main()
