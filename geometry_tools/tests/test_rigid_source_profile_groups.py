import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.rigid_source_profile_groups import rigid_source_profile_groups

class RigidSourceGroupsTest(unittest.TestCase):
    def tube(self,aliases=False):
        theta=np.arange(6)*2*np.pi/6
        points=np.array([[.001*np.cos(t),.002*j,.001*np.sin(t)]for j in range(8)for t in theta]);faces=[]
        for j in range(7):
            for k in range(6):
                a=j*6+k;b=j*6+(k+1)%6;c=a+6;d=b+6;faces.extend([[a,b,c],[b,d,c]])
        faces=np.array(faces)
        if aliases:
            points=np.r_[points,points[19:20]];faces[faces==19]=len(points)-1;faces[25,1]=19
        return Mesh(points,faces,'test:metres')

    def test_sections_are_complete_and_the_two_open_ends_are_fixed(self):
        mesh=self.tube();r=rigid_source_profile_groups(mesh)
        self.assertGreater(len(r['groups']),2)
        for j in range(8):self.assertEqual(len(set(r['group_of_vertex'][j*6:(j+1)*6])),1)
        self.assertTrue(set(r['group_of_vertex'][np.r_[0:6,42:48]])<=set(r['boundary_groups']))
        self.assertEqual(sorted(np.concatenate(r['groups']).tolist()),list(range(len(mesh.vertices))))

    def test_exact_position_aliases_never_acquire_different_owners(self):
        mesh=self.tube(aliases=True);r=rigid_source_profile_groups(mesh)
        self.assertEqual(r['group_of_vertex'][19],r['group_of_vertex'][-1])
        self.assertFalse(r['anatomicalIdentityInferred'])
        self.assertTrue(r['sourceFaceOwnershipComplete'])

    def test_local_core_patches_do_not_bind_a_long_unrecovered_chart(self):
        v=np.array([[0.,0.,0.],[.03,0.,0.],[0.,.03,0.],[0.,0.,.03]])
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        mesh=Mesh(v,f,'test:metres')
        whole=rigid_source_profile_groups(mesh)
        local=rigid_source_profile_groups(mesh,maximum_core_patch_m=.0015)
        self.assertEqual(len(whole['groups']),1)
        self.assertEqual(len(local['groups']),4)
        self.assertTrue(local['sourceFaceOwnershipComplete'])

    def test_local_core_mode_preserves_all_actual_tube_sections(self):
        mesh=self.tube(aliases=True);r=rigid_source_profile_groups(mesh,maximum_core_patch_m=.0015)
        for j in range(8):self.assertEqual(len(set(r['group_of_vertex'][j*6:(j+1)*6])),1)
        self.assertEqual(r['group_of_vertex'][19],r['group_of_vertex'][-1])

if __name__=='__main__':unittest.main()
