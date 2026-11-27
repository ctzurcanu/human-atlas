import itertools,unittest
import numpy as np
from geometry_tools.algorithms.integer_arrangement import IntegerArrangement
from geometry_tools.model import Mesh
from geometry_tools.audit import topology


class IntegerArrangementTests(unittest.TestCase):
    def partition(self,points,tets,order,maximum=False):
        scores=np.eye(len(points))*2-1
        pairs=list(itertools.combinations(range(len(points)),2))
        fields=np.c_[-np.ones(len(points)),np.column_stack([scores[:,i]-scores[:,j]for i,j in pairs])]
        def owner(v):
            return str(int(np.argmax(np.r_[0.,-v[1:len(points)]])))
        a=IntegerArrangement(points,fields,owner)
        p,surfaces,cells=(a.build_maximum(np.array(tets)[order],len(points))if maximum else a.build(np.array(tets)[order]))
        volume=0.
        for name,uses in cells.items():
            f=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
            r=topology(Mesh(p,f,'fixture:metres'))
            self.assertTrue(r['closed_oriented_manifold'],(name,r))
            self.assertGreater(r['oriented_triangle_integral_m3'],0)
            volume+=r['oriented_triangle_integral_m3']
        return a,volume

    def test_dependent_comparison_planes_share_closed_cells(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
        a,v=self.partition(p,[[0,1,2,3]],[0])
        self.assertAlmostEqual(v,1/6)
        self.assertEqual(len(a.affine_index),len(a.points))

    def test_adjacent_tetrahedra_cancel_their_internal_faces_in_either_order(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
        for order in [[0,1],[1,0]]:
            _,v=self.partition(p,[[0,1,2,3],[0,2,1,4]],order)
            self.assertAlmostEqual(v,1/3)

    def test_noninteger_fields_are_rejected(self):
        with self.assertRaises(ValueError):IntegerArrangement(np.eye(3),np.array([[0.],[.5],[1.]]),lambda v:'x')

    def test_direct_maximum_keeps_closed_shared_cells(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
        for order in [[0,1],[1,0]]:
            _,v=self.partition(p,[[0,1,2,3],[0,2,1,4]],order,maximum=True)
            self.assertAlmostEqual(v,1/3)

    def test_generic_material_level_removes_a_four_sector_zero_edge_pinch(self):
        p=np.array([[0.,0,-1],[0,0,1],[1,0,0],[0,1,0],[-1,0,0],[0,-1,0]])
        t=np.array([[0,1,2,3],[0,1,3,4],[0,1,4,5],[0,1,5,2]])
        u=np.array([0,0,1,-1,1,-1])*8+1
        a=IntegerArrangement(p,np.c_[-np.ones(6),u],lambda v:None)
        points,surfaces,cells=a.build_maximum(t,2)
        self.assertTrue(np.array_equal(points[:len(p)],p))
        volume=0.
        for uses in cells.values():
            f=np.concatenate([surfaces[q['surface']]if q['side']==1 else surfaces[q['surface']][:,[0,2,1]]for q in uses])
            result=topology(Mesh(points,f,'fixture:metres'))
            self.assertTrue(result['closed_oriented_manifold'])
            volume+=result['oriented_triangle_integral_m3']
        self.assertAlmostEqual(volume,4/3)

    def test_dominance_pruning_keeps_a_partial_stage_interface_closed(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
        scores=(np.eye(5)*2-1)*128+np.arange(5)[None]
        pairs=list(itertools.combinations(range(5),2))
        fields=np.c_[-np.ones(5),np.column_stack([scores[:,i]-scores[:,j]for i,j in pairs])]
        a=IntegerArrangement(p,fields,lambda v:None)
        a.build_maximum(np.array([[0,1,2,3]]),5)
        points,surfaces,cells=a.build_maximum(np.array([[0,2,1,4]]),5)
        volume=0.
        for uses in cells.values():
            f=np.concatenate([surfaces[q['surface']]if q['side']==1 else surfaces[q['surface']][:,[0,2,1]]for q in uses])
            result=topology(Mesh(points,f,'fixture:metres'))
            self.assertTrue(result['closed_oriented_manifold'])
            volume+=result['oriented_triangle_integral_m3']
        self.assertAlmostEqual(volume,1/3)

    def test_priority_coat_has_its_own_scale_and_keeps_adjacent_tets_closed(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
        t=np.array([[0,1,2,3],[0,2,1,4]])
        # The two planes cross inside the tissue. A priority face must share
        # its actual subdivision with both underlying materials.
        fields=np.c_[-np.ones(5),(4*p[:,1]-1)*1000000,4*p[:,0]-1]
        for order in [[0,1],[1,0]]:
            a=IntegerArrangement(p,fields,lambda v:None)
            a.build_maximum(t[order[:1]],2,priority_field=2,priority_label='fascia')
            points,surfaces,cells=a.build_maximum(t[order[1:]],2,priority_field=2,priority_label='fascia')
            volumes={}
            for name,uses in cells.items():
                faces=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
                r=topology(Mesh(points,faces,'fixture:metres'))
                self.assertTrue(r['closed_oriented_manifold'],(name,r))
                volumes[name]=r['oriented_triangle_integral_m3']
            self.assertAlmostEqual(volumes['fascia'],9/64)
            self.assertAlmostEqual(volumes['region:0'],19/192)
            self.assertAlmostEqual(volumes['region:1'],3/32)
            self.assertAlmostEqual(sum(volumes.values()),1/3)
            self.assertTrue(np.array_equal(points[:len(p)],p))
            for use in cells['fascia']:
                if use.get('neighbor')is not None:continue
                faces=surfaces[use['surface']]
                self.assertFalse(np.any(np.all(points[faces,0]==.25,axis=1)),
                                 'The internal priority boundary must be reciprocal, never a free exterior')

    def test_priority_on_an_existing_tet_face_has_no_zero_volume_tissue(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
        t=np.array([[0,1,2,3],[0,2,1,4]])
        fields=np.c_[-np.ones(5),2*p[:,0]-1,p[:,2]]
        a=IntegerArrangement(p,fields,lambda v:None)
        points,surfaces,cells=a.build_maximum(t,2,priority_field=2,priority_label='fascia')
        volumes={}
        for name,uses in cells.items():
            faces=np.concatenate([surfaces[u['surface']]if u['side']==1 else surfaces[u['surface']][:,[0,2,1]]for u in uses])
            r=topology(Mesh(points,faces,'fixture:metres'))
            self.assertTrue(r['closed_oriented_manifold'],(name,r))
            self.assertGreater(r['oriented_triangle_integral_m3'],0)
            volumes[name]=r['oriented_triangle_integral_m3']
        self.assertAlmostEqual(volumes['fascia'],1/6)
        self.assertAlmostEqual(sum(volumes.values()),1/3)

    def test_priority_cannot_reuse_a_material_comparison_or_label(self):
        p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
        a=IntegerArrangement(p,np.c_[-np.ones(4),np.arange(4),np.arange(4)-1],lambda v:None)
        with self.assertRaises(ValueError):a.build_maximum([[0,1,2,3]],2,priority_field=1)
        with self.assertRaises(ValueError):a.build_maximum([[0,1,2,3]],2,priority_field=2,priority_label='region:0')
