import unittest
from collections import defaultdict
import numpy as np
from geometry_tools.algorithms.shared_triangle_edges import boundary,reconcile_boundaries,refine_surfaces,bisect_surfaces
from geometry_tools.algorithms.finite_tissue_boolean import CommonPositionBank,make_solid
from geometry_tools.audit import topology
from geometry_tools.model import Mesh

class SharedEdgeTests(unittest.TestCase):
    def test_refinement_preserves_area_and_reciprocal_edge(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.5,0.,0.]])
        expected=np.array([[0,1,2]]);actual=np.array([[0,3,2],[3,1,2]])
        plans=defaultdict(set);self.assertEqual(reconcile_boundaries(v,actual,expected,plans),0)
        bank=CommonPositionBank(v);surfaces,parents,changed=refine_surfaces(bank,{'a':expected,'b':expected[:,[0,2,1]]},plans)
        self.assertEqual(changed,2);self.assertEqual(boundary(surfaces['a']),boundary(actual));self.assertEqual(boundary(surfaces['b'][:,[0,2,1]]),boundary(actual))
        xyz=bank.positions[surfaces['a']];self.assertAlmostEqual(np.linalg.norm(np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]),axis=1).sum()/2,.5)
        self.assertTrue(np.array_equal(bank.positions[:len(v)],v));self.assertTrue(np.all(parents['a']==0))
    def test_different_boundary_rejected(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.5,.0001,0.]])
        with self.assertRaises(ValueError):reconcile_boundaries(v,np.array([[0,3,2],[3,1,2]]),np.array([[0,1,2]]),defaultdict(set))
    def test_same_curve_different_orientation_rejected(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        with self.assertRaises(ValueError):reconcile_boundaries(v,np.array([[0,2,1]]),np.array([[0,1,2]]),defaultdict(set))
    def test_all_midpoint_split_patterns_preserve_orientation_area_and_shared_boundary(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[.5,0.,0.],[.5,.5,0.],[0.,.5,0.]])
        face=np.array([[0,1,2]]);edges=[(0,1),(1,2),(0,2)]
        for mask in range(1,8):
            plans={e:3+i for i,e in enumerate(edges)if mask&(1<<i)};surfaces,parents,count=bisect_surfaces({'a':face,'b':face[:,[0,2,1]]},plans);xyz=v[surfaces['a']];normal=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);self.assertTrue(np.all(normal[:,2]>0));self.assertAlmostEqual(normal[:,2].sum()/2,.5);self.assertEqual(boundary(surfaces['a']),boundary(surfaces['b'][:,[0,2,1]]));self.assertEqual(count,2);self.assertTrue(np.all(parents['a']==0))
    def test_global_propagation_preserves_two_closed_finite_cells(self):
        # Two tetrahedra share one tissue chart. A cut inserts an edge control
        # on just one owner's version; all incident exterior charts need it.
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.],[0.,0.,-1.],[.5,0.,0.]])
        interface=np.array([[0,2,1]])
        split=np.array([[0,2,5],[5,2,1]])
        upper=np.array([[0,1,3],[0,3,2],[1,2,3]])
        lower=np.array([[0,4,1],[0,2,4],[1,4,2]])
        plans=defaultdict(set)
        self.assertEqual(reconcile_boundaries(v,interface,split,plans),0.)
        bank=CommonPositionBank(v)
        charts,parents,changed=refine_surfaces(bank,dict(interface=interface,upper=upper,lower=lower),plans)
        self.assertEqual(changed,3)
        self.assertTrue(np.array_equal(bank.positions[:len(v)],v))
        for name,side in [('upper',1),('lower',-1)]:
            shared=charts['interface']if side==1 else charts['interface'][:,[0,2,1]]
            fs=np.concatenate([shared,charts[name]])
            mesh=Mesh(bank.positions,fs,'analytic-shared-tetrahedra')
            check=topology(mesh)
            self.assertTrue(check['closed_oriented_manifold'])
            self.assertEqual(len(check['components']),1)
            self.assertEqual(check['components'][0]['closed_orientable_genus'],0)
            solid,_=make_solid(bank.positions,fs)
            self.assertAlmostEqual(solid.volume(),1/6,places=14)
    def test_thin_parent_does_not_snap_a_new_centroid_to_its_boundary(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,3e-13,0.],[1/3,0.,0.]])
        bank=CommonPositionBank(v)
        charts,parents,_=refine_surfaces(bank,{'thin':np.array([[0,1,2]])},{(0,1):{3}})
        self.assertTrue(np.array_equal(bank.positions,v))
        self.assertEqual(len(charts['thin']),2)
        xyz=bank.positions[charts['thin']];areas=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])[:,2]
        self.assertTrue(np.all(areas>0))
        np.testing.assert_allclose(areas.sum(),3e-13,rtol=4e-15,atol=0)
        self.assertEqual(boundary(charts['thin']),boundary(np.array([[0,3,2],[3,1,2]])))

if __name__=='__main__':unittest.main()
