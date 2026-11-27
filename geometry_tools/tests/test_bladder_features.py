"""A shared ridge moves wall thickness without filling a real opening."""
import unittest
import numpy as np

from geometry_tools.bladder_features import back_graph, TrigonalRelief
from geometry_tools.model import Mesh
from geometry_tools.harmonic_wall import mesher_coordinates


class BladderFeatureTests(unittest.TestCase):
    def test_mesher_frame_retains_all_source_bits_at_anatomical_scales(self):
        rng=np.random.default_rng(57)
        p=rng.uniform([-.04,.82,-.05],[.04,1.3,.13],(4096,3))
        local,origin,scale=mesher_coordinates(p)
        np.testing.assert_array_equal(local/scale+origin,p)
        self.assertEqual(scale,1024.)
        self.assertEqual(origin.tolist(),[0.,1.,0.])

    def fixture(self):
        xy = np.array([[-.04,.88],[.04,.88],[.04,.95],[-.04,.95]])
        v = np.r_[np.c_[xy,np.full(4,-.024)],np.c_[xy,np.full(4,-.020)]]
        outer = np.array([[0,2,1],[0,3,2]])
        inner = np.array([[4,5,6],[4,6,7]])
        mesh = Mesh(v, np.r_[outer,inner], 'test:metres')
        a = np.arange(16)*2*np.pi/16
        rims = {name:np.c_[x+.002*np.cos(a),y+.002*np.sin(a),np.full(16,-.02)]
                for name,x,y in [('left',.021,.925),('right',-.021,.925),('neck',0,.89)]}
        return TrigonalRelief(mesh, inner, outer, rims), rims

    def test_missing_sheet_remains_a_miss_instead_of_extrapolation(self):
        v = np.array([[0,0,-.002],[1,0,-.002],[0,1,-.002]])
        z, parents = back_graph(v, np.array([[0,1,2]]), np.array([[.2,.2],[.8,.8]]))
        self.assertEqual(z[0], -.002)
        self.assertTrue(np.isnan(z[1]))
        self.assertEqual(parents.tolist(), [0,-1])

    def test_one_map_thickens_the_wall_and_preserves_outer_sheet_and_ostia(self):
        relief,rims = self.fixture()
        p = np.array([[0,.925,-.024],[0,.925,-.022],[0,.925,-.020],
                      [0,.8965,-.020],[0,.925,.02]])
        delta,bound = relief.displacement(p)
        np.testing.assert_allclose(delta, [0,.0005,.001,.0008,0], atol=2e-17)
        self.assertTrue(np.all(bound > 0))
        for ring in rims.values():
            np.testing.assert_array_equal(relief.displacement(ring)[0], 0)

    def test_relief_is_monotone_through_tissue_and_the_fluid_falloff(self):
        relief,_ = self.fixture()
        z = np.linspace(-.03,.0,5001)
        p = np.c_[np.zeros(len(z)),np.full(len(z),.925),z]
        moved = z+relief.displacement(p)[0]
        self.assertTrue(np.all(np.diff(moved) > 0))
        self.assertGreater(float(np.min(np.diff(moved)/np.diff(z))), .62)


if __name__ == '__main__':
    unittest.main()
