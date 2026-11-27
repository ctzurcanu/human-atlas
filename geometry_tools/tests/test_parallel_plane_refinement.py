import unittest
import numpy as np
from geometry_tools.algorithms.parallel_plane_refinement import refine_parallel_planes


class ParallelPlaneRefinementTest(unittest.TestCase):
    def test_closed_tetrahedron_retains_exact_surface_and_source_incidence(self):
        v=np.array([[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        r=refine_parallel_planes(v,f,[.2,.4,.7])
        self.assertTrue(np.array_equal(r['vertices'][:len(v)],v))
        replay=(v[r['vertex_source_edges']]*r['vertex_source_weights'][:,:,None]).sum(1)
        self.assertLess(np.max(abs(replay-r['vertices'])),1e-14)
        e=np.sort(r['faces'][:,[0,1,1,2,2,0]].reshape(-1,2),axis=1)
        _,counts=np.unique(e,axis=0,return_counts=True)
        self.assertTrue(np.all(counts==2));self.assertEqual(len(r['vertices'])-len(counts)+len(r['faces']),2)
        self.assertEqual(set(r['parent_faces']),set(range(len(f))))

    def test_plane_at_existing_vertex_adds_no_false_cap_or_duplicate_triangle(self):
        v=np.array([[0.,0.,0.],[1.,1.,0.],[0.,2.,0.],[-1.,1.,0.]])
        f=np.array([[0,1,2],[0,2,3]])
        r=refine_parallel_planes(v,f,[1.])
        self.assertEqual(len(r['faces']),4)
        self.assertTrue(np.all(r['vertices'][:,2]==0))
        self.assertEqual(len(np.unique(np.sort(r['faces'],axis=1),axis=0)),len(r['faces']))

    def test_clockwise_faces_keep_their_orientation(self):
        v=np.array([[0.,0.,0.],[0.,2.,0.],[1.,1.,0.]])
        r=refine_parallel_planes(v,[[0,1,2]],[.5,1.,1.5])
        q=r['vertices'][r['faces']]
        self.assertTrue(np.all(np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])[:,2]<0))


if __name__=='__main__':unittest.main()
