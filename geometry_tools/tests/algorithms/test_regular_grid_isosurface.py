import unittest
import numpy as np
from geometry_tools.algorithms.regular_grid_isosurface import extract
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

class RegularGridSurfaceTests(unittest.TestCase):
    def test_sphere_is_one_positive_closed_volume(self):
        q=np.linspace(-1.4,1.4,23);x,y,z=np.meshgrid(q,q,q,indexing='ij');d=np.sqrt(x*x+y*y+z*z)-1
        v,f,_=extract(d,[q[0]]*3,q[1]-q[0]);t=topology(Mesh(v,f,'test:metres'))
        self.assertTrue(t['closed_oriented_manifold']);self.assertEqual([x['closed_orientable_genus']for x in t['components']],[0]);self.assertGreater(t['oriented_triangle_integral_m3'],3.9)
    def test_through_lumen_stays_open(self):
        q=np.linspace(-1.4,1.4,25);x,y,z=np.meshgrid(q,q,q,indexing='ij');d=np.maximum(np.sqrt(x*x+y*y+z*z)-1,.35-np.sqrt(x*x+y*y))
        v,f,_=extract(d,[q[0]]*3,q[1]-q[0]);t=topology(Mesh(v,f,'test:metres'))
        self.assertTrue(t['closed_oriented_manifold']);self.assertEqual([x['closed_orientable_genus']for x in t['components']],[1]);self.assertTrue(np.all(np.linalg.norm(v[:,:2],axis=1)>.30))
    def test_explicit_cell_budget_cannot_expand(self):
        q=np.linspace(-1.4,1.4,12);x,y,z=np.meshgrid(q,q,q,indexing='ij')
        with self.assertRaisesRegex(ValueError,'budget'):extract(np.sqrt(x*x+y*y+z*z)-1,[q[0]]*3,q[1]-q[0],maximum_active_cells=1)

if __name__=='__main__':unittest.main()
