import unittest
import numpy as np
from geometry_tools.algorithms.selective_plane_refinement import refine_selected_plane


class SelectivePlaneRefinementTests(unittest.TestCase):
    def test_local_square_cut_propagates_only_shared_edge(self):
        p=np.array([[-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]])
        f=np.array([[0,1,2],[0,2,3]])
        r=refine_selected_plane(p,f,normal=[1,0,0],coordinate=0.,active_faces=[True,False])
        self.assertEqual(len(r['positions']),6);self.assertEqual(len(r['faces']),5)
        self.assertEqual(r['evidence']['propagated_faces'],1)
        q=r['positions'][r['faces']]
        replay=np.einsum('tki,tij->tkj',r['corner_barycentrics'],p[f[r['parent_faces']]])
        np.testing.assert_array_equal(q,replay)
        self.assertAlmostEqual(float(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])[:,2].sum()/2),4.)
        edge=np.sort(np.concatenate([r['faces'][:,[0,1]],r['faces'][:,[1,2]],r['faces'][:,[2,0]]]),axis=1)
        unique,count=np.unique(edge,axis=0,return_counts=True)
        for a,b in unique[count==1]:
            self.assertTrue(np.any(np.all(np.abs(r['positions'][[a,b],:2])==1,axis=0)))

    def test_one_tetrahedral_face_cut_keeps_entire_boundary_closed(self):
        p=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
        r=refine_selected_plane(p,f,normal=[1,0,0],coordinate=.37,active_faces=[1,0,0,0])
        e=np.sort(np.concatenate([r['faces'][:,[0,1]],r['faces'][:,[1,2]],r['faces'][:,[2,0]]]),axis=1)
        _,counts=np.unique(e,axis=0,return_counts=True);self.assertTrue(np.all(counts==2))
        q=r['positions'][r['faces']]
        volume=np.einsum('ni,ni->n',q[:,0],np.cross(q[:,1],q[:,2])).sum()/6
        self.assertAlmostEqual(float(volume),1/6,places=15)
        np.testing.assert_array_equal(r['positions'][:len(p)],p)

    def test_oblique_plane_and_near_corner_classification(self):
        p=np.array([[-1.,-1.,0.],[1.,-1.,0.],[1.,1.,0.],[-1.,1.,0.]])
        f=np.array([[0,1,2],[0,2,3]])
        r=refine_selected_plane(p,f,normal=[1,-1,0],coordinate=.3)
        q=r['positions'][r['faces']]
        self.assertTrue(np.all(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])[:,2]>0))
        d=q[:,:,0]-q[:,:,1]-.3
        self.assertTrue(np.all((d.min(1)>=-1e-14)|(d.max(1)<=1e-14)))
        p=np.array([[1e-16,1.,0.],[-1.,0.,0.],[1.,0.,0.]])
        r=refine_selected_plane(p,[[0,1,2]],normal=[1,0,0],coordinate=0.,plane_tolerance_m=1e-12)
        self.assertEqual(len(r['faces']),2);np.testing.assert_array_equal(r['positions'][:3],p)

    def test_empty_active_region_is_literal(self):
        p=np.array([[-1.,0.,0.],[1.,0.,0.],[0.,1.,0.]])
        r=refine_selected_plane(p,[[0,1,2]],normal=[1,0,0],coordinate=0.,active_faces=[False])
        np.testing.assert_array_equal(r['positions'],p);np.testing.assert_array_equal(r['faces'],[[0,1,2]])

if __name__=='__main__':unittest.main()
