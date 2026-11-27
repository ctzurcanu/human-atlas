import unittest
import numpy as np
from geometry_tools.algorithms.closed_curve_correspondence import closed_curve_correspondence,align_existing_subdivision
from geometry_tools.model import Mesh
from geometry_tools.source_port_dependencies import check_source_mouth
from geometry_tools.algorithms.source_mouth_index import extract_mouths,match_mouths,declared_seam_slab


class ClosedCurveTests(unittest.TestCase):
    def setUp(self):
        self.p=np.array([[0.,0.,0.],[.001,0,0],[.001,.001,0],[0,.001,0]])
        self.q=np.insert(self.p,1,[.0003,0,0],axis=0)

    def test_unequal_subdivisions_of_the_same_curve(self):
        r=closed_curve_correspondence(self.p,np.roll(self.q[::-1],2,axis=0))
        self.assertLess(r['maximumGap_m'],1e-12)
        self.assertTrue(r['wholeCurveBound'])

    def test_deviation_between_coarse_corners_is_included(self):
        q=self.q.copy();q[1,2]=.0001
        r=closed_curve_correspondence(self.p,q[::-1])
        self.assertGreaterEqual(r['maximumGap_m'],.0001)

    def test_wrong_order_and_direction_cannot_be_scrambled(self):
        r=closed_curve_correspondence(self.p,self.q)
        self.assertGreater(r['maximumGap_m'],.0005)

    def test_complete_mesh_boundaries_support_unequal_counts(self):
        a=Mesh(self.p,np.array([[0,1,2],[0,2,3]]),'fixture')
        b=Mesh(self.q,np.array([[0,3,1],[1,3,2],[0,4,3]]),'fixture')
        r=check_source_mouth(a,b,range(4),range(5))
        self.assertTrue(r['passed']);self.assertTrue(r['curveCorrespondence']['wholeCurveBound'])
        mouths=extract_mouths(a,owner='a',system='venous')['mouths']+extract_mouths(b,owner='b',system='venous')['mouths']
        index=match_mouths(mouths)
        self.assertEqual(len(index['passed']),1)
        self.assertTrue(index['passed'][0]['curveCorrespondence']['wholeCurveBound'])
        with self.assertRaises(ValueError):check_source_mouth(a,b,range(4),[0,2,3,4])

    def test_degenerate_and_nonfinite_curves_rejected(self):
        with self.assertRaises(ValueError):closed_curve_correspondence(self.p,self.p[[0,0,2]])
        q=self.q.copy();q[1,0]=np.nan
        with self.assertRaises(ValueError):closed_curve_correspondence(self.p,q)

    def test_alignment_keeps_fine_corners_and_moves_only_the_declared_rims(self):
        fine=self.q.copy();fine[1,2]=2e-6
        coarse=self.p[::-1]+[0,0,1e-6]
        a,b,proof=align_existing_subdivision(fine,coarse,maximum_motion_m=3e-6)
        np.testing.assert_array_equal(a[[0,2,3,4]],fine[[0,2,3,4]])
        self.assertLess(closed_curve_correspondence(a,b)['maximumGap_m'],1e-12)
        self.assertEqual(len(a),5);self.assertEqual(len(b),4)
        with self.assertRaises(ValueError):align_existing_subdivision(fine,coarse,maximum_motion_m=1e-6)

    def test_nearest_points_cannot_replace_cyclic_corner_order(self):
        with self.assertRaises(ValueError):align_existing_subdivision(self.q,self.p,maximum_motion_m=1e-6)

    def test_unequal_source_seam_certificate_has_a_strict_rounding_bound(self):
        def band(rim,height,reverse=False):
            n=len(rim);v=np.r_[rim,rim+[0,0,height]];faces=[]
            for i in range(n):
                j=(i+1)%n;faces.extend([[i,j,n+i],[j,n+j,n+i]])
            f=np.asarray(faces)
            return Mesh(v,f[:,[0,2,1]]if reverse else f,'fixture')
        a,b=band(self.p,.001),band(self.q,-.001,True)
        proof=declared_seam_slab(a,b,range(4),range(5),0,0,maximum_width_m=1e-8)
        self.assertTrue(proof['qualified']);self.assertFalse(proof['anatomicalAcceptance'])
        shifted=Mesh(b.vertices+[0,0,1e-7],b.faces,b.frame)
        self.assertFalse(declared_seam_slab(a,shifted,range(4),range(5),0,0,maximum_width_m=1e-8)['qualified'])
        curved=self.p.copy();curved[2,2]=.0001
        refined=np.insert(curved,1,[.0003,0,0],axis=0)
        first,second=band(curved,.001),band(refined,-.001,True)
        proof=declared_seam_slab(first,second,range(4),range(5),0,0,maximum_width_m=1e-8)
        self.assertTrue(proof['qualified']);self.assertEqual(proof['method'],'curved-rim-separating-axis')


if __name__=='__main__':unittest.main()
