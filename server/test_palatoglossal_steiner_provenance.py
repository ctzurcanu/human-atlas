"""Actual five-face interior witness and strict negative boundary fixtures."""
import sys,unittest
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT));sys.path.insert(0,str(Path(__file__).parent))
from palatoglossal_original_provenance import interior_steiner_control
from geometry_tools.algorithms.exact_source_intersections import contains
from palatoglossal_incidence_fixture import qvec
W=ROOT/'.local-models/geometry-workbench/tongue-completion'
class ActualInteriorControl(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        z=np.load(W/'palatoglossal-direct-arch-steiner-witness-v1/actual-star.npz');cls.p=z['point'];cls.star=z['star_triangles'];cls.parent=z['source_triangle']
    def test_actual_five_face_star_exactly_preserves_binary_encoding(self):
        self.assertEqual(len(self.star),5);q=interior_steiner_control(self.p,self.star,[(5,405)],self.parent)
        self.assertEqual(tuple(float(x)for x in q),tuple(self.p));self.assertTrue(contains(q,tuple(qvec(p)for p in self.parent)))
    def test_actual_boundary_star_does_not_qualify(self):
        with self.assertRaisesRegex(ValueError,'not a closed oriented link'):
            interior_steiner_control(self.p,self.star[:-1],[(5,405)],self.parent)
    def test_actual_four_face_star_uses_exact_rounding_cell(self):
        z=np.load(W/'palatoglossal-direct-arch-steiner-witness-v2/actual-star.npz');q=interior_steiner_control(z['point'],z['star_triangles'],[(5,405)],z['source_triangle'])
        self.assertEqual(tuple(float(x)for x in q),tuple(z['point']));self.assertTrue(contains(q,tuple(qvec(p)for p in z['source_triangle'])))
    def test_source_plane_outside_rounding_cell_is_rejected(self):
        point=np.array([0.,0.,1.]);rim=np.array([[1.,0.,1.],[0.,1.,1.],[-1.,0.,1.],[0.,-1.,1.]]);star=np.array([[point,rim[i],rim[(i+1)%4]]for i in range(4)])
        with self.assertRaisesRegex(ValueError,'misses the identical Float64 rounding cell'):
            interior_steiner_control(point,star,[(1,0)],np.array([[-2.,-2.,0.],[2.,-2.,0.],[0.,2.,0.]]))
    def test_multiple_original_planes_remain_rejected(self):
        with self.assertRaisesRegex(ValueError,'multiple original parents'):
            interior_steiner_control(self.p,self.star,[(5,405),(5,406)],self.parent)
    def test_parent_boundary_is_not_interior_steiner(self):
        point=np.array([0.,0.,0.]);rim=np.array([[1.,0.,0.],[0.,1.,0.],[-1.,0.,0.],[0.,-1.,0.]]);star=np.array([[point,rim[i],rim[(i+1)%4]]for i in range(4)])
        with self.assertRaisesRegex(ValueError,'on an original parent boundary'):
            interior_steiner_control(point,star,[(1,0)],np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]]))
    def test_distinct_disconnected_stars_do_not_qualify(self):
        with self.assertRaises(ValueError):interior_steiner_control(self.p,np.r_[self.star,self.star],[(5,405)],self.parent)
    def test_actual_seven_face_disconnected_star_is_rejected(self):
        z=np.load(W/'palatoglossal-direct-arch-steiner-witness-v3/actual-star.npz')
        with self.assertRaisesRegex(ValueError,'disconnected link cycles'):
            interior_steiner_control(z['point'],z['star_triangles'],[(5,405)],z['source_triangle'])
if __name__=='__main__':unittest.main()
