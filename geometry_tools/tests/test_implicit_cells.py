import unittest
import numpy as np
from geometry_tools.implicit_cells import grid_tetrahedra,Arrangement,has_occupied_volume
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.distance import triangle_pairs,exact_transverse_crossing


class CommonFieldTests(unittest.TestCase):
    def test_contact_polygons_have_no_tissue_volume_but_thin_cells_are_retained(self):
        planar=np.array([[0.,0,0],[1,0,0],[1,1,0],[0,1,0]])
        self.assertFalse(has_occupied_volume(planar))
        self.assertTrue(has_occupied_volume(np.r_[planar,[[.5,.5,1e-17]]]))
        self.assertFalse(has_occupied_volume(np.array([[0.,0,0],[1,0,0],[0,1,0]])))
        arrangement=Arrangement(planar,np.zeros((4,1)),lambda x:'contact')
        arrangement.emit([[0,1,2,3],[0,3,2],[0,2,1]],'contact')
        self.assertEqual(arrangement.boundaries,{})
        self.assertEqual(arrangement.zero_volume_leaves,1)

    def test_finite_shell_and_regional_cells_are_reciprocal(self):
        p,t,_=grid_tetrahedra([-1.4]*3,[1.4]*3,.36);r=np.linalg.norm(p,axis=1)
        v=np.array([r-1.1,r-.7,p[:,0]-.12]).T
        a=Arrangement(p,v,lambda x:None if x[0]>=0 else 'shell'if x[1]>=0 else'left'if x[2]<0 else'right',vertex_snap_m=1e-12)
        p,s,c=a.build(t)
        self.assertEqual(set(c),{'left','right','shell'})
        for name,uses in c.items():
            f=np.concatenate([s[q['surface']]if q['side']==1 else s[q['surface']][:,[0,2,1]]for q in uses])
            q=audit(Mesh(p,f,'fixture'),require_closed=True,expected_genus=[0]if name!='shell'else[0,0])
            self.assertTrue(q['passed'],q['failures'])
            for use in uses:
                if use['neighbor']:
                    peer=next(x for x in c[use['neighbor']]if x['surface']==use['surface'])
                    self.assertEqual(peer['side'],-use['side']);self.assertEqual(peer['neighbor'],name)

    def test_exact_predicate_rejects_near_coplanar_shared_edge_false_hit(self):
        a=np.array([[.066686631318,1.220721101536,-.020947392635],[.067090037721,1.220105629798,-.021429364764],[.066857537464,1.220326529178,-.021041948167]])
        b=np.array([a[0],[.066857496991,1.220326815907,-.021041948167],a[1]])
        self.assertFalse(exact_transverse_crossing(a,b))
        self.assertEqual(triangle_pairs(a[None],b[None])['relation'][0],'touching')
        a=np.array([[-1.,0,0],[1,0,0],[0,-1,0]]);b=np.array([[0,-.5,-.01],[0,.5,-.01],[0,0,1.]])
        self.assertTrue(exact_transverse_crossing(a,b));self.assertEqual(triangle_pairs(a[None],b[None])['relation'][0],'crossing')

    def test_sparse_grid_keeps_identical_material_extent(self):
        lower,upper=[-1.4]*3,[1.4]*3
        p,t,_=grid_tetrahedra(lower,upper,.4);q,u,_=grid_tetrahedra(lower,upper,.4,domain_distance=lambda p:np.linalg.norm(p,axis=1)-.8)
        material=lambda x:None if x[0]>=0 else'one'
        def build(v,f):
            a=Arrangement(v,(np.linalg.norm(v,axis=1)-.8)[:,None],material);v,s,c=a.build(f)
            return Mesh(v,next(iter(s.values())),'fixture')
        self.assertEqual(build(p,t).surface_sha,build(q,u).surface_sha)


if __name__=='__main__':unittest.main()
