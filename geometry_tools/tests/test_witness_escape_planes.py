import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.witness_escape_planes import witnessed_face_planes


class WitnessEscapePlaneTests(unittest.TestCase):
    def test_interior_triangle_crossing_missed_by_vertices_midpoints_and_centroid_gets_a_complete_plane(self):
        v=np.array([[0.,0,0],[.01,0,0],[0,.01,0]]);f=np.array([[0,1,2]])
        box=np.array([[x,y,z]for x in [.001,.0015]for y in [.002,.0025]for z in [-.0001,.0001]])
        from scipy.spatial import ConvexHull
        hull=ConvexHull(box);obstacle=Mesh(box,hull.simplices,'test:metres')
        samples=np.r_[v,(v+np.roll(v,-1,axis=0))*.5,v.mean(0,keepdims=True)]
        self.assertFalse(np.any(np.all((samples>=box.min(0))&(samples<=box.max(0)),axis=1)))
        hits,proof=pair_contacts(Mesh(v,f,'test:metres'),obstacle,0,max_tests=1000,max_pairs=1000)
        self.assertTrue(proof['complete_within_pair_budget']);self.assertGreater(len(hits['pairs']),0)
        guards=witnessed_face_planes(v+[0,0,.001],v,f,obstacle,hits['pairs'],part_id='bone')
        self.assertEqual(len(guards),1);g=guards[0];n=np.array(g['normal']);point=np.array(g['point_m'])
        shift=max(0,g['clearance_m']-float(((v-point)@n).min()))*n
        self.assertTrue(np.all((v+shift-point)@n>=g['clearance_m']-1e-15))
        hits,_=pair_contacts(Mesh(v+shift,f,'test:metres'),obstacle,0,max_tests=1000,max_pairs=1000)
        self.assertEqual(len(hits['pairs']),0)
        reverse=witnessed_face_planes(v+[0,0,.001],v,f,obstacle,np.array([[0,0]]),part_id='outer-skin',desired_side='away-from-reference')
        np.testing.assert_allclose(reverse[0]['normal'],[0.,0.,-1.])
        with self.assertRaises(ValueError):witnessed_face_planes(v,v,f,obstacle,np.array([[0,0]]),part_id='bone')


if __name__=='__main__':unittest.main()
