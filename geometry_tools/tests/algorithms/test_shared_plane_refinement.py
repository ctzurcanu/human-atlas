import unittest
import numpy as np
from geometry_tools.algorithms.shared_plane_refinement import refine_at_axis_plane,safe_plane_coordinate,lift_refined_positions


class SharedPlaneRefinementTests(unittest.TestCase):
    def test_shared_edge_is_split_once_with_exact_corner_ancestry(self):
        p=np.array([[-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]])
        f=np.array([[0,1,2],[0,2,3]])
        r=refine_at_axis_plane(p,f,axis=0,coordinate=0.)
        q=r['positions'][r['faces']]
        replay=np.einsum('tki,tij->tkj',r['corner_barycentrics'],p[f[r['parent_faces']]])
        np.testing.assert_array_equal(q,replay)
        self.assertEqual(len(r['positions']),7)
        self.assertEqual(len(r['faces']),6)
        self.assertEqual(int(np.all(r['positions']==[0.,0.,0.],axis=1).sum()),1)
        area=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])[:,2]/2
        self.assertTrue(np.all(area>0));self.assertAlmostEqual(float(area.sum()),4.)
        self.assertTrue(np.all((q[:,:,0].max(1)<=0)|(q[:,:,0].min(1)>=0)))

    def test_plane_through_one_corner_preserves_orientation(self):
        p=np.array([[-1.,0.,0.],[0.,1.,0.],[1.,0.,0.]])
        r=refine_at_axis_plane(p,[[0,1,2]],axis=0,coordinate=0.)
        q=r['positions'][r['faces']]
        self.assertEqual(len(q),2)
        self.assertTrue(np.all(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])[:,2]<0))
        self.assertEqual(len(r['positions']),4)

    def test_tetrahedron_stays_closed_with_equal_oriented_volume(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
        r=refine_at_axis_plane(p,f,axis=0,coordinate=.37)
        edges=np.sort(np.concatenate([r['faces'][:,[0,1]],r['faces'][:,[1,2]],r['faces'][:,[2,0]]]),axis=1)
        _,counts=np.unique(edges,axis=0,return_counts=True);self.assertTrue(np.all(counts==2))
        q=r['positions'][r['faces']]
        volume=np.einsum('ni,ni->n',q[:,0],np.cross(q[:,1],q[:,2])).sum()/6
        self.assertAlmostEqual(float(volume),1/6,places=15)

    def test_unsafe_near_vertex_plane_is_rejected(self):
        with self.assertRaises(ValueError):
            refine_at_axis_plane([[0,0,0],[1,0,0],[0,1,0]],[[0,1,2]],axis=0,coordinate=1e-12)

    def test_sparse_interval_and_outside_planes_remain_valid(self):
        self.assertEqual(safe_plane_coordinate([0.,.5],.35,search_radius_m=.01),.35)
        self.assertEqual(safe_plane_coordinate([0.,.5],.8,search_radius_m=.01),.8)

    def test_unsafe_plane_without_available_gap_is_rejected(self):
        with self.assertRaises(ValueError):
            safe_plane_coordinate([0.,1.],1e-12,search_radius_m=1e-8)

    def test_oblique_frame_subdivision_retains_literal_world_controls(self):
        p=np.array([[-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]])
        f=np.array([[0,1,2],[0,2,3]])
        angle=.2;c,s=np.cos(angle),np.sin(angle)
        matrix=np.array([[c,-s,0.],[s,c,0.],[0.,0.,1.]])
        world=p@matrix.T+[.021,1.24,.01]
        r=refine_at_axis_plane(p,f,axis=0,coordinate=0.)
        lifted=lift_refined_positions(world,f,r)
        np.testing.assert_array_equal(lifted[:len(p)],world)
        np.testing.assert_allclose(lifted,r['positions']@matrix.T+[.021,1.24,.01],atol=3e-16)

if __name__=='__main__':unittest.main()
