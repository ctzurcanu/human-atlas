"""Regression evidence for native footprint clipping and degenerate cutters."""
import sys,unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'scripts'))
from ta98_projected_clip import clip_triangle,clip_affine,subtract_convex,signed_area,clean_polygon
class ProjectedClip(unittest.TestCase):
 def test_inner_cutter_has_disjoint_area(self):
  p=np.array([[0.,0.],[1,0],[1,1],[0,1]]);q=np.array([[.2,.2],[.8,.2],[.8,.8],[.2,.8]])
  pieces=subtract_convex(p,q);self.assertAlmostEqual(sum(abs(signed_area(x)) for x in pieces),.64,places=12)
  for i,a in enumerate(pieces):
   for b in pieces[i+1:]:
    for j in range(1,len(b)-1):self.assertLess(abs(signed_area(clip_triangle(a,b[[0,j,j+1]]))),1e-12)
 def test_zero_area_cutter_does_not_duplicate_surface(self):
  p=np.array([[0.,0.],[1,0],[1,1],[0,1]]);q=np.array([[.5,0],[.5,.5],[.5,1]])
  self.assertAlmostEqual(sum(abs(signed_area(x)) for x in subtract_convex(p,q)),1,places=12)
 def test_tiny_duplicate_edge_has_no_extra_area(self):
  p=np.array([[0.,0.],[.01,0],[.01,.01],[0,.01]]);q=np.array([[.002,.002],[.008,.002],[.008,.002+1e-18],[.008,.008],[.002,.008]])
  self.assertAlmostEqual(sum(abs(signed_area(x)) for x in subtract_convex(p,q)),.000064,places=14)
 def test_affine_clipping_invariant_under_positive_coefficient_scale(self):
  p=np.array([[0.,0.],[1,0],[1,1],[0,1]])
  self.assertAlmostEqual(abs(signed_area(clip_affine(p,[1,0,-.5]))),abs(signed_area(clip_affine(p,[1e-12,0,-.5e-12]))),places=12)
 def test_native_scale_arithmetic_duplicates_do_not_create_cutter_edges(self):
  p=np.array([[.010789279825985435,-.02756073512136935],
              [.010789279825985433,-.027560735121369372],
              [.010789279825985432,-.027560735121369362],
              [.011643877252936361,-.027347087860107425],
              [.01088971000951744,-.027202055697911474]])
  q=clean_polygon(p)
  self.assertEqual(len(q),3)
  self.assertLess(abs(signed_area(q)-signed_area(p)),1e-18)
  self.assertTrue(np.array_equal(q[1:],p[-2:]))
if __name__=='__main__':unittest.main()
