import unittest
from unittest.mock import patch
from types import SimpleNamespace
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.monotone_context_clearance import monotone_context_clearance
from geometry_tools.tests.test_sparse_group_clearance import cube


class MonotoneContextTests(unittest.TestCase):
    def setUp(self):
        self.bone=cube([0.,0.,0.],[.01,.01,.01])
        self.tree=BVH(self.bone);self.index=SurfaceIndex(self.bone)
        self.mesh=Mesh(np.array([[.0099,.004,.004],[.0101,.006,.004],[.0101,.005,.006]]),
                       np.array([[0,1,2]]),'test:metres')

    def test_clearance_keeps_complete_rigid_section(self):
        q,proof=monotone_context_clearance(self.mesh,np.arange(3),[np.arange(3)],[],self.tree,self.index)
        self.assertTrue(proof['complete'],proof)
        self.assertGreater(proof['initialTargetPairs'],0)
        self.assertEqual(proof['finalTargetPairs'],0)
        np.testing.assert_allclose(q-self.mesh.vertices,np.broadcast_to(q[0]-self.mesh.vertices[0],q.shape),atol=1e-15)
        self.assertTrue(all(s['after']<s['before'] for s in proof['acceptedSteps']))
        self.assertEqual(proof['finalSelfPairs'],0)

    def test_protected_group_cannot_be_moved(self):
        q,proof=monotone_context_clearance(self.mesh,np.arange(3),[np.arange(3)],[0],self.tree,self.index)
        self.assertFalse(proof['complete'])
        np.testing.assert_array_equal(q,self.mesh.vertices)
        self.assertEqual(proof['initialTargetPairs'],proof['finalTargetPairs'])

    def test_open_skin_sheet_moves_a_crossing_section_inward(self):
        # The x-positive face is an open skin patch, not a closed solid.
        # An occupancy query on this sheet cannot establish its inside.
        faces=self.bone.faces[np.all(self.bone.triangles[:,:,0]==.01,axis=1)]
        skin=Mesh(self.bone.vertices,faces,self.bone.frame)
        q,proof=monotone_context_clearance(self.mesh,np.arange(3),[np.arange(3)],[],
              BVH(skin),SurfaceIndex(skin),target_side='interior-sheet')
        self.assertTrue(proof['complete'],proof)
        self.assertLess(float(q[:,0].max()),.01)
        self.assertEqual(proof['finalTargetPairs'],0)
        self.assertFalse(proof['openSkinContainmentEstablished'])
        self.assertTrue(proof['sourceFaceOrientationGuard'])

    def test_unknown_target_side_is_rejected(self):
        with self.assertRaises(ValueError):
            monotone_context_clearance(self.mesh,np.arange(3),[],[],self.tree,self.index,
                                       target_side='infer-from-overlap')

    def test_geodesic_spread_does_not_move_a_nearby_disconnected_branch(self):
        faces=self.bone.faces[np.all(self.bone.triangles[:,:,0]==.01,axis=1)]
        skin=Mesh(self.bone.vertices,faces,self.bone.frame)
        other=self.mesh.vertices+[-.001,0,.0001]
        mesh=Mesh(np.r_[self.mesh.vertices,other],np.array([[0,1,2],[3,4,5]]),self.mesh.frame)
        q,proof=monotone_context_clearance(mesh,np.arange(6),[np.arange(3),np.arange(3,6)],[],
              BVH(skin),SurfaceIndex(skin),target_side='interior-sheet',geodesic_spread_m=.01,
              corotate_sections=True)
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_array_equal(q[3:],other)
        self.assertLess(q[:3,0].max(),.01)

    def test_corotation_requires_a_declared_field(self):
        with self.assertRaises(ValueError):
            monotone_context_clearance(self.mesh,np.arange(3),[],[],self.tree,self.index,
                                       corotate_sections=True)

    def test_skin_without_crossings_is_not_clear_when_entirely_outside(self):
        faces=self.bone.faces[np.all(self.bone.triangles[:,:,0]==.01,axis=1)]
        skin=Mesh(self.bone.vertices,faces,self.bone.frame)
        mesh=Mesh(self.mesh.vertices+[.002,0,0],self.mesh.faces,self.mesh.frame)
        q,proof=monotone_context_clearance(mesh,np.arange(3),[np.arange(3)],[],
              BVH(skin),SurfaceIndex(skin),target_side='interior-sheet')
        self.assertEqual(proof['finalTargetPairs'],0)
        self.assertFalse(proof['complete'])
        self.assertGreater(proof['finalSampledSkinViolation']['maximumViolation_m'],.001)
        np.testing.assert_array_equal(q,mesh.vertices)

    def test_skin_side_descent_rejects_an_undeclared_skin_target(self):
        with self.assertRaises(ValueError):
            monotone_context_clearance(self.mesh,np.arange(3),[],[],self.tree,self.index,
                                       skin_penetration_descent=True)

    def test_skin_side_descent_has_explicit_decreasing_energy(self):
        faces=self.bone.faces[np.all(self.bone.triangles[:,:,0]==.01,axis=1)]
        skin=Mesh(self.bone.vertices,faces,self.bone.frame)
        q,proof=monotone_context_clearance(self.mesh,np.arange(3),[np.arange(3)],[],
              BVH(skin),SurfaceIndex(skin),target_side='interior-sheet',skin_penetration_descent=True)
        self.assertTrue(proof['complete'],proof)
        self.assertLess(proof['finalSampledSkinViolation']['squaredViolation_m2'],
                        proof['initialSampledSkinViolation']['squaredViolation_m2'])
        for step in proof['acceptedSteps']:
            self.assertLessEqual(step['after'],step['before'])
            self.assertLess(step['sampledSkinViolationAfter']['squaredViolation_m2'],
                            step['sampledSkinViolationBefore']['squaredViolation_m2'])

    def test_penetration_can_decrease_before_the_crossing_count_changes(self):
        faces=self.bone.faces[np.all(self.bone.triangles[:,:,0]==.01,axis=1)]
        skin=Mesh(self.bone.vertices,faces,self.bone.frame)
        mesh=Mesh(np.array([[.009,.004,.004],[.0118,.006,.004],[.0118,.005,.006]]),
                  self.mesh.faces,self.mesh.frame)
        unchanged,strict=monotone_context_clearance(mesh,np.arange(3),[np.arange(3)],[],
              BVH(skin),SurfaceIndex(skin),target_side='interior-sheet',maximum_motion_m=.001)
        np.testing.assert_array_equal(unchanged,mesh.vertices)
        moved,descent=monotone_context_clearance(mesh,np.arange(3),[np.arange(3)],[],
              BVH(skin),SurfaceIndex(skin),target_side='interior-sheet',maximum_motion_m=.001,
              skin_penetration_descent=True)
        self.assertGreater(len(descent['acceptedSteps']),0,descent)
        self.assertTrue(any(s['before']==s['after']for s in descent['acceptedSteps']))
        self.assertLess(descent['finalSampledSkinViolation']['maximumViolation_m'],
                        descent['initialSampledSkinViolation']['maximumViolation_m'])
        self.assertFalse(descent['complete'],'The remaining outside wall must stay unaccepted')
        self.assertLessEqual(np.linalg.norm(moved-mesh.vertices,axis=1).max(),.001+1e-15)

    def test_connected_patch_clears_adjacent_penetrating_controls_together(self):
        mesh=Mesh(np.array([[.0098,.004,.004],[.0101,.006,.004],
                            [.0101,.006,.006],[.0098,.004,.006]]),
                  np.array([[0,1,2],[0,2,3]]),'test:metres')
        q,proof=monotone_context_clearance(mesh,np.arange(4),[],[],self.tree,self.index,
                                          connected_contact_groups=True)
        self.assertTrue(proof['complete'],proof)
        self.assertEqual(proof['finalTargetPairs'],0)
        self.assertEqual(proof['finalSelfPairs'],0)
        np.testing.assert_allclose(q-mesh.vertices,np.broadcast_to(q[0]-mesh.vertices[0],q.shape),atol=1e-15)
        self.assertTrue(proof['connectedContactGroups'])

    def test_late_final_check_cannot_report_complete(self):
        mesh=Mesh(self.mesh.vertices+[.003,0.,0.],self.mesh.faces,self.mesh.frame)
        ticks=iter([0.,13.])
        with patch('geometry_tools.algorithms.monotone_context_clearance.time',
                   SimpleNamespace(monotonic=lambda:next(ticks))):
            q,proof=monotone_context_clearance(mesh,np.arange(3),[np.arange(3)],[],self.tree,self.index)
        self.assertEqual(proof['finalTargetPairs'],0)
        self.assertFalse(proof['complete'])
        self.assertTrue(proof['timedOut'])
        np.testing.assert_array_equal(q,mesh.vertices)


if __name__=='__main__':unittest.main()
