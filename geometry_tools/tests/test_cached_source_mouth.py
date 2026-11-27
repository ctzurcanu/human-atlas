import unittest

import numpy as np

from geometry_tools.algorithms.cached_source_mouth import SourceMouthChecks
from geometry_tools.model import Mesh
from geometry_tools.source_port_dependencies import check_source_mouth
from geometry_tools.tests.test_declared_mouth_batch import DeclaredMouthBatchTests


class CachedMouthTests(unittest.TestCase):
    def test_literal_proof_and_one_index_per_immutable_mesh(self):
        meshes,pair=DeclaredMouthBatchTests().fixtures(gap=5e-7)
        a,b=meshes.values();first,second=pair['firstRing'],pair['secondRing'];cache=SourceMouthChecks()
        expected=check_source_mouth(a,b,first,second)
        for _ in range(4):self.assertEqual(cache.check(a,b,first,second),expected)
        self.assertEqual(cache.geometry_builds,2)

    def test_frame_and_incomplete_curve_errors_match(self):
        meshes,pair=DeclaredMouthBatchTests().fixtures();a,b=meshes.values();cache=SourceMouthChecks()
        for first,second in [(pair['firstRing'][:-1],pair['secondRing']),([6,6,8,9,10,11],pair['secondRing'])]:
            for fn in [check_source_mouth,cache.check]:
                with self.assertRaisesRegex(ValueError,'complete external source mouth'):fn(a,b,first,second)
        other=Mesh(b.vertices,b.faces,'other:metres',b.attributes,b.source)
        for fn in [check_source_mouth,cache.check]:
            with self.assertRaisesRegex(ValueError,'same explicit coordinate frame'):fn(a,other,pair['firstRing'],pair['secondRing'])

    def test_new_mesh_has_independent_topology_and_tolerance(self):
        meshes,pair=DeclaredMouthBatchTests().fixtures(gap=5e-6);a,b=meshes.values();cache=SourceMouthChecks()
        self.assertFalse(cache.check(a,b,pair['firstRing'],pair['secondRing'])['passed'])
        self.assertTrue(cache.check(a,b,pair['firstRing'],pair['secondRing'],tolerance_m=1e-5)['passed'])
        p=b.vertices.copy();p[0]=p[1];changed=Mesh(p,b.faces,b.frame,b.attributes,b.source)
        for fn in [check_source_mouth,cache.check]:
            with self.assertRaises(ValueError):fn(a,changed,pair['firstRing'],pair['secondRing'])

    def test_complete_unequal_boundary_curve_proof_is_literal(self):
        meshes,pair=DeclaredMouthBatchTests().fixtures(gap=5e-7);a,b=meshes.values()
        p=np.vstack([b.vertices,(b.vertices[0]+b.vertices[1])*.5])
        faces=np.vstack([[[0,12,6],[12,1,6]],b.faces[1:]])
        refined=Mesh(p,faces,b.frame)
        ring=[0,12,1,2,3,4,5]
        expected=check_source_mouth(a,refined,pair['firstRing'],ring)
        self.assertTrue(expected['passed']);self.assertIn('curveCorrespondence',expected)
        self.assertEqual(SourceMouthChecks().check(a,refined,pair['firstRing'],ring),expected)


if __name__=='__main__':unittest.main()
