"""Portable reservoir fixtures: openings remain fluid passages, never solid parts."""
from pathlib import Path
import tempfile
import unittest

import manifold3d as md
import numpy as np

from geometry_tools.bladder_carriers import recover
from geometry_tools.model import Mesh


class BladderCarrierTests(unittest.TestCase):
    def test_three_port_reservoir_replays_outer_minus_fluid(self):
        outer=md.Manifold.sphere(1.,16);inner=md.Manifold.sphere(.8,16)
        port=md.Manifold.cylinder(2.,.15,.15,8)
        fluid=inner+port+port.rotate([0,90,0])+port.rotate([-90,0,0])
        wall=outer-fluid;q=wall.to_mesh64()
        source=Mesh(np.asarray(q.vert_properties)[:,:3],np.asarray(q.tri_verts,int),'fixture')
        with tempfile.TemporaryDirectory()as td:
            out=Path(td)/'carriers';r=recover(source,out)
            self.assertEqual(r['mouth_count'],3);self.assertTrue(r['source_boundary_chain_exact'])
            self.assertFalse(r['computational_disks_for_renderer']);self.assertFalse(r['atlas_modified'])
            self.assertTrue(all(q['passed']for q in r['checks'].values()))
            self.assertAlmostEqual(r['volumes_m3']['wall'],wall.volume(),places=12)
            self.assertAlmostEqual(r['volumes_m3']['outer']-r['volumes_m3']['fluid'],wall.volume(),places=12)


if __name__=='__main__':unittest.main()
