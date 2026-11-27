"""The graph's safety inflation must not silently replace curve constraints."""
import tempfile
from pathlib import Path
import unittest
import numpy as np
from geometry_tools.algorithms.communicating_vascular_network import InteriorRouter

class CurveClearanceMasks(unittest.TestCase):
    def test_separate_curve_mask_retains_default_and_forbidden_endpoints(self):
        vertices=np.array([[x,y,z]for x in [-.01,.01]for y in [-.01,.01]for z in [-.01,.01]])
        faces=np.array([[0,1,3],[0,3,2],[4,6,7],[4,7,5],[0,4,5],[0,5,1],[2,3,7],[2,7,6],[0,2,6],[0,6,4],[1,5,7],[1,7,3]])
        index=np.array([[x,y,z]for x in range(-8,9)for y in range(-4,5)for z in range(-4,5)])
        points=index*.001
        with tempfile.TemporaryDirectory()as directory:
            cache=Path(directory)/'grid.npz'
            np.savez(cache,index=index,points=points,distance=np.full(len(points),.005))
            router=InteriorRouter.from_cache(vertices[faces],cache,pitch=.001)
            graph_mask=np.abs(points[:,0])>.0065
            a,b=np.array([-.0075,.0013,.0017]),np.array([.0075,.0013,.0017])
            with self.assertRaisesRegex(ValueError,'constrained clearance'):
                router.route(a,b,forbidden=graph_mask)
            curve=router.route(a,b,forbidden=graph_mask,curve_forbidden=np.zeros(len(points),bool))
            np.testing.assert_array_equal(curve[[0,-1]],[a,b])
            self.assertTrue(np.all(np.abs(curve)<.01))
            with self.assertRaisesRegex(ValueError,'constrained clearance'):
                router.route(a,b,forbidden=graph_mask,curve_forbidden=graph_mask)

if __name__=='__main__':unittest.main()
