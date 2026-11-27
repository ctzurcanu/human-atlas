import unittest
import numpy as np
from geometry_tools.algorithms.scalar_surface_clip import clip_scalar_surface,interpolate_clip_controls


def area(v,f):
    p=v[f];return np.linalg.norm(np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]),axis=1).sum()/2


class ScalarSurfaceClipTests(unittest.TestCase):
    def test_analytic_triangle_partition_and_no_cap(self):
        v=np.array([[0,0,0],[2,0,0],[0,2,0]],float);r=clip_scalar_surface(v,np.array([[0,1,2]]),v[:,0],1)
        self.assertAlmostEqual(area(r['positions'],r['positive_faces']),.5)
        self.assertAlmostEqual(area(r['positions'],r['negative_faces']),1.5)
        np.testing.assert_array_equal(r['positions'][:3],v);self.assertEqual(r['caps_added'],0)
        self.assertEqual(set(map(tuple,r['positions'][3:])),{(1.,0.,0.),(1.,1.,0.)})

    def test_neighbors_share_one_exact_cut_control(self):
        v=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0]],float);f=np.array([[0,1,2],[0,2,3]])
        r=clip_scalar_surface(v,f,v[:,0],.5)
        self.assertEqual(int(np.all(r['positions']==[.5,.5,0],axis=1).sum()),1)
        self.assertEqual(len(r['new_edge_parents']),3)
        self.assertAlmostEqual(area(r['positions'],r['positive_faces']),.5)
        self.assertAlmostEqual(area(r['positions'],r['negative_faces']),.5)

    def test_corresponding_wall_interface_uses_the_same_subdivision(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0]],float);r=clip_scalar_surface(v,np.array([[0,1,2]]),v[:,0],.25)
        inner=v+[0,0,-.1];q=interpolate_clip_controls(inner,r['new_edge_parents'],r['fractions'])
        np.testing.assert_allclose(q,r['positions']+[0,0,-.1],atol=0,rtol=0)

    def test_original_isovalue_vertex_is_not_duplicated(self):
        v=np.array([[0,0,0],[1,0,0],[0,1,0]],float);r=clip_scalar_surface(v,np.array([[0,1,2]]),np.array([0,1,-1.]),0)
        self.assertEqual(len(r['positions']),4)
        self.assertAlmostEqual(area(r['positions'],r['positive_faces']),.25)
        self.assertAlmostEqual(area(r['positions'],r['negative_faces']),.25)

    def test_unresolved_face_and_bad_inputs_are_rejected(self):
        v=np.eye(3);f=np.array([[0,1,2]])
        for s,level in [(np.zeros(3),0),(np.zeros(2),1),(np.array([0,np.nan,1]),.5)]:
            with self.assertRaises(ValueError):clip_scalar_surface(v,f,s,level)
        with self.assertRaises(ValueError):clip_scalar_surface(v,f.astype(float),np.arange(3),.5)


if __name__=='__main__':unittest.main()
