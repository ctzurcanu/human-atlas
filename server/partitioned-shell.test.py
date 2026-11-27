"""Analytic finite-volume tests; never reads or modifies anatomical models."""
import sys
import unittest
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT/'scripts'))
from geometry_tools.algorithms.partitioned_shell import build_partitioned_shell, validate_reciprocal_ownership
from ta98_geometry_checks import topology
from ta98_triangle_intersections import sweep_hits


def fixture():
    # Octahedron halves; explicitly specified equatorial tissue interface.
    outer = np.array([[.01,0,0],[0,.01,0],[-.01,0,0],[0,-.01,0],
                      [0,0,.01],[0,0,-.01]], dtype='<f4')
    inner = outer*np.float32(.8)
    faces = np.array([[0,1,4],[1,2,4],[2,3,4],[3,0,4],
                      [1,0,5],[2,1,5],[3,2,5],[0,3,5]], int)
    return outer, inner, faces, np.array([0]*4+[1]*4), ('top','bottom')


class FinitePartitionedShell(unittest.TestCase):
    def assembly(self):
        return build_partitioned_shell(*fixture(), interior_interfaces=[
            {'a':'top','b':'bottom','faces':[[6,8,7],[6,9,8]]}])

    def test_closed_distinct_volumes_tile_carrier_without_duplicate_interfaces(self):
        assembly = self.assembly()
        self.assertTrue(validate_reciprocal_ownership(assembly))
        total = 0
        for cell in assembly.cells:
            triangles = assembly.cell_faces(cell)
            result = topology(assembly.positions, triangles, tolerance=0)
            for key in ['openEdges','nonmanifoldEdges','inconsistentWindingEdges',
                        'duplicateFaces','degenerateFaces']:
                self.assertEqual(result[key], 0, (cell,key))
            self.assertEqual(result['invalidCycleOrPathVertexLinks'], [])
            self.assertEqual(len(result['components']), 1)
            self.assertGreater(result['signedVolumeMetersCubed'], 0)
            hits = sweep_hits(assembly.positions[triangles], assembly.positions[triangles], True)
            self.assertEqual(hits['properOrCoplanarIntersections'], 0)
            total += result['signedVolumeMetersCubed']
        outer,inner,faces,_,_ = fixture()
        self.assertAlmostEqual(total, topology(outer,faces,tolerance=0)['signedVolumeMetersCubed'], places=18)
        self.assertTrue(np.array_equal(assembly.positions[:len(outer)], outer))
        self.assertTrue(np.array_equal(assembly.positions[len(outer):], inner))
        # Four external/inner patches, one rim, one interior interface.
        self.assertEqual(len(assembly.surfaces), 6)
        exterior = assembly.exterior_faces()
        self.assertEqual(sorted(map(tuple, exterior)), sorted(map(tuple, faces)))

    def test_missing_internal_interface_remains_open_and_cannot_pass_topology(self):
        assembly = build_partitioned_shell(*fixture())
        self.assertTrue(validate_reciprocal_ownership(assembly))
        self.assertEqual(topology(assembly.positions,assembly.cell_faces('core:top'))['openEdges'],4)

    def test_input_arrays_remain_exact_and_output_bank_is_readonly(self):
        args = fixture()
        before = [value.tobytes() for value in args[:4]]
        assembly = build_partitioned_shell(*args)
        self.assertEqual(before,[value.tobytes() for value in args[:4]])
        with self.assertRaises(ValueError):
            assembly.positions[0] = 0

    def test_rejects_invalid_ownership_open_whole_and_exterior_internal_caps(self):
        outer,inner,faces,owners,regions = fixture()
        with self.assertRaises(ValueError):
            build_partitioned_shell(outer,inner,faces,owners[:-1],regions)
        with self.assertRaises(ValueError):
            build_partitioned_shell(outer,inner,faces[:-1],owners[:-1],regions)
        with self.assertRaises(ValueError):
            build_partitioned_shell(*fixture(),interior_interfaces=[
                {'a':'top','b':'bottom','faces':[[0,1,2]]}])


if __name__ == '__main__':
    unittest.main()
