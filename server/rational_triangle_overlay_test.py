from fractions import Fraction as F
from pathlib import Path
import sys,unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from geometry_tools.algorithms.rational_triangle_overlay import RationalSurfaceBank,source_weights


class ExactOverlayTests(unittest.TestCase):
    def setUp(self):
        self.bank=RationalSurfaceBank([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        for i in range(3):self.bank.register(i,{i:F(1)},maximum_replay_error_m=0.)
        self.parent=(0,1,2)

    def test_crossed_subdivisions_preserve_area_and_point_identity(self):
        b=self.bank;a=b.add({0:F(1,2),1:F(1,2)});c=b.add({0:F(2,3),2:F(1,3)})
        first=[[0,a,2],[a,1,2]];second=[[0,1,c],[c,1,2]];area=F(0);shared=set()
        for x in first:
            for y in second:
                p=b.triangle_intersection(x,y,self.parent);r=b.triangle_intersection(y,x,self.parent)
                self.assertEqual(set(p),set(r))
                if not p:continue
                area+=b.signed_area2(p,self.parent)
                for f in b.triangulate(p,self.parent):self.assertGreater(b.signed_area2(f,self.parent),0)
                shared.update(p)
        self.assertEqual(area,F(1))
        crossing=b.add({0:F(2,5),1:F(2,5),2:F(1,5)})
        self.assertIn(crossing,shared)

    def test_touching_has_no_positive_area(self):
        a=self.bank.add({0:F(1,2),1:F(1,2)})
        self.assertEqual(self.bank.triangle_intersection([0,a,2],[a,1,2],self.parent),[])

    def test_exact_edge_identity_and_unrelated_edge_rejection(self):
        a=self.bank.add({0:F(1,3),1:F(2,3)})
        self.assertEqual(self.bank.segment_fraction(a,0,1),F(2,3))
        self.assertIsNone(self.bank.segment_fraction(a,0,2))

    def test_declared_simplex_and_replay_bound_are_enforced(self):
        w=source_weights([.5,0.,0.],[0,1],np.array(self.bank.rows))
        self.assertEqual(w,{0:F(1,2),1:F(1,2)})
        with self.assertRaisesRegex(ValueError,'arithmetic bound'):
            bank=RationalSurfaceBank([[0.,0.,.1],[1.,0.,0.],[0.,1.,0.],[.25,.25,.9]])
            bank.register(3,{0:F(1,2),1:F(1,4),2:F(1,4)},maximum_replay_error_m=.01)


if __name__=='__main__':unittest.main()
