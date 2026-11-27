import unittest

import numpy as np

from geometry_tools.algorithms.declared_mouth_batch import reconcile_declared_mouths
from geometry_tools.model import Mesh


class DeclaredMouthBatchTests(unittest.TestCase):
    def fixtures(self, gap=5e-6):
        t=np.arange(6)*np.pi/3
        radial=np.c_[np.cos(t),np.sin(t),np.zeros(6)]
        faces=[]
        for i in range(6):
            j=(i+1)%6
            faces.extend([[i,j,i+6],[j,j+6,i+6]])
        meshes={}
        for name,z in [('first',[0,.001]),('second',[.001+gap,.002])]:
            p=np.concatenate([radial*.0002+[0,0,k] for k in z])
            meshes[name]=Mesh(p,np.array(faces),'test:metres',dict(
                normals=np.rint(np.tile(radial,(2,1))*32767).astype('<i2'),
                normalSigns=np.ones(12,dtype='i1')),{})
        return meshes,dict(first='first',second='second',firstRing=list(range(6,12)),secondRing=list(range(6)))

    def solve(self, meshes, pairs):
        return reconcile_declared_mouths(meshes,pairs,maximum_gap_m=3e-5,
            maximum_motion_m=2e-5,maximum_relative_area_change=.025,
            maximum_facet_rotation_degrees=10)

    def test_complete_common_rim_no_caps_or_distal_movement(self):
        meshes,pair=self.fixtures()
        result,proof=self.solve(meshes,[pair])
        np.testing.assert_array_equal(result['first'].vertices[6:],result['second'].vertices[:6])
        np.testing.assert_array_equal(result['first'].vertices[:6],meshes['first'].vertices[:6])
        np.testing.assert_array_equal(result['second'].vertices[6:],meshes['second'].vertices[6:])
        for name in meshes:
            np.testing.assert_array_equal(result[name].faces,meshes[name].faces)
            np.testing.assert_array_equal(result[name].attributes['normalSigns'],meshes[name].attributes['normalSigns'])
        self.assertEqual(proof['lumenDisksCreated'],0)

    def test_far_source_ends_are_rejected(self):
        meshes,pair=self.fixtures(gap=.0001)
        with self.assertRaises(ValueError):self.solve(meshes,[pair])

    def test_double_owned_port_is_rejected(self):
        meshes,pair=self.fixtures()
        with self.assertRaisesRegex(ValueError,'overlapping'):
            self.solve(meshes,[pair,pair])


if __name__=='__main__':unittest.main()
