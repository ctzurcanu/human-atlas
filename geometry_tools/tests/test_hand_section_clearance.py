import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.commands.clear_hand_profiles import section_proposal
from geometry_tools.algorithms.profile_core_transport import profile_translation_weights,transport_profile_cores
from geometry_tools.commands.filter_hand_components import encoded_attribute


class HandSectionClearanceTests(unittest.TestCase):
    def test_component_merge_keeps_packed_normals_at_full_precision(self):
        normals=np.array([[1.,0,0],[0,-1,0],[1,1,0]])
        encoded=encoded_attribute('normals',normals,np.dtype('int16'))
        np.testing.assert_array_equal(encoded[:2],[[32767,0,0],[0,-32767,0]])
        self.assertGreater(encoded[2,0],23000)
        np.testing.assert_array_equal(encoded_attribute('normals',encoded,np.dtype('int16')),encoded)

    def test_full_tissue_influence_replays_core_and_section_aliases(self):
        angle=np.arange(6)*np.pi/3;v=np.concatenate([np.c_[np.full(6,x),np.cos(angle),np.sin(angle)]for x in [0,1,2]])
        faces=[]
        for k in range(2):
            for j in range(6):
                a=k*6+j;b=k*6+(j+1)%6;faces.extend([[a,b,a+6],[b,b+6,a+6]])
        # Middle station belongs to the intervening core, not a prescribed ring.
        profiles=[np.arange(6),np.arange(12,18)];f=np.array(faces)
        v=np.r_[v,v[:1]];f[f==0]=18;physical=np.r_[np.arange(18),0]
        weights=profile_translation_weights(v,f,profiles,physical_ids=physical)
        delta=np.array([[.02,.03,.04],[-.01,.02,.03]])
        target=np.array([v[p].mean(0)for p in profiles])+delta
        actual=transport_profile_cores(v,f,profiles,target,np.repeat(np.eye(3)[None],2,axis=0),physical_ids=physical)
        np.testing.assert_allclose(v+weights@delta,actual['vertices'],atol=1e-14)
        np.testing.assert_array_equal(weights[0],weights[18])
        np.testing.assert_allclose(weights.sum(1),1,atol=1e-14)

    def test_joint_course_moves_complete_sections_between_two_surfaces(self):
        bone_vertices=np.array([[-.02,-.01,-.02],[.02,-.01,-.02],[.02,.01,-.02],[-.02,.01,-.02],
            [-.02,-.01,.001],[.02,-.01,.001],[.02,.01,.001],[-.02,.01,.001]])
        box_faces=np.array([[0,2,1],[0,3,2],[4,5,6],[4,6,7],[0,1,5],[0,5,4],
            [1,2,6],[1,6,5],[2,3,7],[2,7,6],[3,0,4],[3,4,7]])
        bone=SurfaceIndex(Mesh(bone_vertices,box_faces,'fixture-metres'))
        skin=SurfaceIndex(Mesh(np.array([[-.03,-.03,.005],[.03,-.03,.005],[.03,.03,.005],[-.03,.03,.005]]),
            np.array([[0,1,2],[0,2,3]]),'fixture-metres'))
        angle=np.arange(6)*np.pi/3
        v=np.concatenate([np.c_[np.full(6,x),.001*np.cos(angle),.0015+.001*np.sin(angle)]for x in [-.01,0,.01]])
        profiles=np.arange(18).reshape(3,6);faces=[]
        for k in range(2):
            for j in range(6):
                a=k*6+j;b=k*6+(j+1)%6;c=a+6;d=b+6
                faces.extend([[a,b,c],[b,d,c]])
        target,proof=section_proposal(v,np.array(faces),profiles,np.array([[0,1],[1,2]]),bone,skin)
        self.assertTrue(proof['complete'],proof)
        translated=v.copy()
        for i,p in enumerate(profiles):translated[p]+=target[i]-v[p].mean(0)
        self.assertGreaterEqual(translated[:,2].min(),.00119999)
        self.assertLessEqual(translated[:,2].max(),.00480001)
        for p in profiles:
            np.testing.assert_allclose(np.linalg.norm(translated[p][:,None]-translated[p][None,:],axis=2),
                np.linalg.norm(v[p][:,None]-v[p][None,:],axis=2),atol=1e-14)


if __name__=='__main__':unittest.main()
