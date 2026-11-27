import unittest
import numpy as np
from geometry_tools.algorithms.coupled_mask_partition import partition_mask
from geometry_tools.algorithms.finite_tissue_boolean import make_solid
from geometry_tools.tests.test_protected_assembly_carve import ProtectedCarveTests


class CoupledMaskTests(unittest.TestCase):
    def test_void_and_owned_region_preserve_the_protected_neighbour(self):
        p, surfaces, cells, normals, uvs, cutter = ProtectedCarveTests().fixture()
        for region in [None, 'ring']:
            with self.subTest(region=region):
                r = partition_mask(p, surfaces, cells, normals, uvs, ['a'],
                    cutter.vert_properties[:, :3], np.asarray(cutter.tri_verts), region=region)
                for name, expected in [('a', .875), ('b', 1.)] + ([('ring', .125)] if region else []):
                    faces = np.concatenate([r['surfaces'][u['surface']] if u['side'] == 1
                        else r['surfaces'][u['surface']][:, [0, 2, 1]] for u in r['cells'][name]])
                    solid, _ = make_solid(r['positions'], faces)
                    self.assertAlmostEqual(solid.volume(), expected, places=12)
                    edges = np.r_[faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]]
                    self.assertTrue(np.all(np.unique(np.sort(edges, axis=1), axis=0, return_counts=True)[1] == 2))
                np.testing.assert_array_equal(r['positions'][:len(p)], p)
                self.assertTrue(r['evidence']['protected_source_regions_preserved'])
                self.assertEqual(r['evidence']['classifier_cap_contributions'], 0)

    def test_virtual_caps_must_not_supply_tissue_faces(self):
        p, surfaces, cells, normals, uvs, cutter = ProtectedCarveTests().fixture()
        with self.assertRaisesRegex(ValueError, 'virtual end cap'):
            partition_mask(p, surfaces, cells, normals, uvs, ['a'],
                cutter.vert_properties[:, :3], np.asarray(cutter.tri_verts),
                cap_faces=set(range(len(cutter.tri_verts))))


if __name__ == '__main__':
    unittest.main()
