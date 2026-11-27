"""Source-boundary recovery where a tissue shares its host exterior."""
import unittest
import numpy as np
from geometry_tools.algorithms.source_face_complement import complement_source_faces,conform_surface_boundary_edges,reduce_encoded_surface_chains


class SourceFaceComplementTests(unittest.TestCase):
    def setUp(self):
        self.v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.5,0.,0.],[0.,.5,0.]])
        self.f=np.array([[0,1,2]])

    @staticmethod
    def area(v,f):
        p=v[f];return float(np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1).sum()/2)

    def test_common_exterior_with_no_contact_is_retained_literally(self):
        v,f,parents,proof=complement_source_faces(self.v,self.v[:3],self.f,np.empty((0,3),int),np.empty(0,int))
        self.assertTrue(np.array_equal(v,self.v));self.assertTrue(np.array_equal(f,self.f))
        self.assertEqual(parents.tolist(),[0]);self.assertEqual(proof['anatomicalCapsAdded'],False)

    def test_near_edge_cut_retains_its_narrow_free_strip_at_declared_precision(self):
        from collections import defaultdict
        from geometry_tools.algorithms.shared_triangle_edges import reconcile_boundaries
        v=np.array([[.04325157031416893,.04839220270514488,.010123497806489468],
                    [.043175000697374344,.04880306124687195,.009721224196255207],
                    [.043224986642599106,.04869589954614639,.009841755963861942],
                    [.04323330225419927,.04860088825188028,.009929898100923791],
                    [.043184164141745876,.04875389188245178,.009769366161747957],
                    [.04323868969825481,.048539353272853226,.009986985118530264]])
        original=np.array([[1,2,0]]);contact=np.array([[2,5,3],[1,2,4],[2,3,4]])
        q,free,_,proof=complement_source_faces(v,v,original,contact,np.zeros(3,int),planar_scale=1e9)
        combined,_=conform_surface_boundary_edges(q,{'patch':np.r_[contact,free]},
            {'one':[dict(surface='patch',side=1)]})
        self.assertLessEqual(reconcile_boundaries(q,combined['patch'],original,defaultdict(set)),2e-12)
        np.testing.assert_array_equal(q[:len(v)],v)
        self.assertAlmostEqual(self.area(q,free)+self.area(q,contact),self.area(v,original),places=18)
        self.assertEqual(proof['localPlanarScale'],1e9)

    def test_partial_contact_preserves_exact_source_area_without_double_cover(self):
        contact=np.array([[0,3,4]])
        v,f,parents,_=complement_source_faces(self.v,self.v[:3],self.f,contact,np.array([0]))
        self.assertAlmostEqual(self.area(v,f)+self.area(v,contact),.5,places=13)
        self.assertTrue(np.array_equal(v,self.v));self.assertTrue(np.all(parents==0))
        self.assertEqual(set(f.ravel()),{1,2,3,4})

    def test_whole_contact_has_no_free_face_and_wrong_parent_plane_fails(self):
        _,f,_,_=complement_source_faces(self.v,self.v[:3],self.f,self.f,np.array([0]))
        self.assertEqual(f.size,0)
        wrong=self.v.copy();wrong[3,2]=.01
        with self.assertRaisesRegex(ValueError,'source parent plane'):
            complement_source_faces(wrong,self.v[:3],self.f,np.array([[0,3,4]]),np.array([0]))

    def test_hanging_edge_refines_both_interface_owners_without_filling_holes(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[.5,0.,0.]])
        surfaces={'interface':np.array([[0,1,2]]),'other':np.array([[0,3,4],[4,3,1],[1,3,2],[2,3,0]])}
        cells={'one':[dict(surface='interface',side=1),dict(surface='other',side=1)],'two':[dict(surface='interface',side=-1)]}
        out,proof=conform_surface_boundary_edges(v,surfaces,cells)
        self.assertEqual(len(out['interface']),2)
        self.assertEqual(proof['positiveAreaCapsAdded'],0);self.assertEqual(proof['newControls'],0)
        # Both uses reference the one canonical chart. The genuinely open
        # reverse owner stays open; the helper must not invent its missing wall.
        self.assertEqual(len(out['other']),4);self.assertEqual(len(out['interface']),2)
        faces=np.r_[out['interface'],out['other']]
        edges=np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]
        _,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
        self.assertTrue(np.all(counts==2))

    def test_open_chart_chain_cancels_only_exact_zero_and_opposed_fragments(self):
        faces=np.array([[0,1,2],[0,1,3],[0,3,4],[0,4,3]])
        out,proof=reduce_encoded_surface_chains(self.v,{'chart':faces})
        self.assertTrue(np.array_equal(out['chart'],self.f))
        self.assertEqual(proof['chart']['exactZeroAreaTrianglesRemoved'],1)
        self.assertEqual(proof['chart']['opposedCoincidentTrianglesCancelled'],2)
        with self.assertRaisesRegex(ValueError,'repeated positive-area'):
            reduce_encoded_surface_chains(self.v,{'chart':np.r_[self.f,self.f]})

    def test_thin_face_uses_the_actual_hanging_edge_without_reversing_a_child(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[.5,1e-12,0.],
                    [0.,0.,-1.],[.25,.5e-12,0.]])
        surfaces={'base':np.array([[0,1,2]]),
                  'sides':np.array([[0,4,3],[4,2,3],[1,0,3],[2,1,3]])}
        cells={'body':[dict(surface='base',side=1),dict(surface='sides',side=1)]}
        out,proof=conform_surface_boundary_edges(v,surfaces,cells)
        self.assertEqual(len(out['base']),2)
        self.assertEqual(proof['newControls'],0)
        p=v[out['base']]
        normal=np.cross(v[1]-v[0],v[2]-v[0])
        self.assertTrue(np.all(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])@normal>0))
        self.assertAlmostEqual(self.area(v,out['base']),self.area(v,surfaces['base']),places=25)
        faces=np.r_[out['base'],out['sides']]
        edges=np.r_[faces[:,[0,1]],faces[:,[1,2]],faces[:,[2,0]]]
        _,counts=np.unique(np.sort(edges,axis=1),axis=0,return_counts=True)
        self.assertTrue(np.all(counts==2))

    def test_nearby_open_patch_does_not_establish_an_opposite_boundary_chain(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],
                    [.5,0.,0.],[.6,-.1,0.],[.4,-.1,0.]])
        surfaces={'first':np.array([[0,1,2]]),'separate':np.array([[3,4,5]])}
        cells={'open':[dict(surface='first',side=1),dict(surface='separate',side=1)]}
        out,proof=conform_surface_boundary_edges(v,surfaces,cells)
        for sid in surfaces:self.assertTrue(np.array_equal(out[sid],surfaces[sid]))
        self.assertEqual(proof['edgeRefinements'],[])

    def test_oblique_source_edge_backtrack_does_not_create_a_reversed_sliver(self):
        v=np.array([[.05021355673670769,1.1761772185564041,-.019676774740219116],
                    [.05164555832743645,1.178322859108448,-.020692821592092514],
                    [.0502481684088707,1.1763610169291496,-.019793877378106117],
                    [.050232889977870676,1.1762798838725703,-.019742185435444395]])
        f=np.array([[0,1,2]]);contact=np.array([[0,1,3]])
        q,free,_,proof=complement_source_faces(v,v[:3],f,contact,np.array([0]))
        self.assertTrue(np.array_equal(q,v))
        self.assertEqual(len(free),1)
        self.assertLess(abs(self.area(q,free)+self.area(q,contact)-self.area(v,f)),1e-16)
        self.assertTrue(proof['numericalCollinearContourVisitsRemoved'])
        original=np.cross(v[1]-v[0],v[2]-v[0]);p=q[free]
        self.assertTrue(np.all(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0])@original>0))


if __name__=='__main__':unittest.main()
