"""Interleaved edge sampling closes by subdivision, never by a cap."""
import unittest
import numpy as np
from geometry_tools.algorithms.source_face_complement import conform_surface_boundary_edges
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

class InterleavedBoundary(unittest.TestCase):
    def test_closed_tetrahedron_with_two_different_edge_samples(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[.3,0.,0.],[.7,0.,0.]])
        # The same A-B edge is sampled at P on one facet and Q on the other.
        f=np.array([[0,2,4],[4,2,1],[0,5,3],[5,1,3],[0,3,2],[1,2,3]])
        cells={'body':[dict(surface='outer',side=1,neighbor=None)]}
        before=topology(Mesh(v,f,'test:metres'));self.assertEqual(before['open_edges'],4)
        surfaces,proof=conform_surface_boundary_edges(v,{'outer':f},cells,allow_collinear_cycles=True)
        after=topology(Mesh(v,surfaces['outer'],'test:metres'))
        self.assertTrue(after['closed_oriented_manifold']);self.assertAlmostEqual(after['oriented_triangle_integral_m3'],1/6)
        self.assertEqual(len(surfaces['outer']),8);self.assertEqual(proof['positiveAreaCapsAdded'],0)
        self.assertEqual(proof['newControls'],0);self.assertFalse(proof['coordinatesMoved'])
        self.assertEqual(len(proof['collinearBoundaryCycleProofs']),1)

    def test_an_actual_opening_remains_open(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        f=np.array([[0,1,3],[0,3,2],[1,2,3]])
        cells={'body':[dict(surface='outer',side=1,neighbor=None)]}
        surfaces,proof=conform_surface_boundary_edges(v,{'outer':f},cells,allow_collinear_cycles=True)
        self.assertTrue(np.array_equal(surfaces['outer'],f));self.assertEqual(proof['collinearBoundaryCycleProofs'],[])
        self.assertEqual(topology(Mesh(v,surfaces['outer'],'test:metres'))['open_edges'],3)

if __name__=='__main__':unittest.main()
