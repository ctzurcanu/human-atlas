import unittest

import numpy as np

from geometry_tools.algorithms.source_control_clearance import source_control_basis
from geometry_tools.algorithms.source_control_clearance import clear_source_controls
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.tests.test_sparse_group_clearance import cube


class SourceControlBasisTests(unittest.TestCase):
    def test_actual_crossing_face_guard_moves_the_entire_rigid_triangle_to_one_side(self):
        v=np.array([[.014,.004,.004],[.016,.006,.004],[.015,.005,.006]])
        sheet=SurfaceIndex(Mesh(np.array([[.015,-.01,-.01],[.015,.02,-.01],[.015,.02,.02],[.015,-.01,.02]]),
                                np.array([[0,1,2],[0,2,3]]),'atlas-world'))
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
            SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
            recipient_surfaces=[(sheet,False,'unsigned-sheet','fixed sheet')],
            recipient_plane_guards=[dict(partId='fixed sheet',sourceRaw=[0,1,2],point_m=[.015,0,0],normal=[1.,0.,0.])],
            rigid_translation=True)
        self.assertTrue(proof['complete'],proof);self.assertGreater(q[:,0].min(),.01519)
        self.assertGreaterEqual(proof['minimumActualRecipientGap_m'],.0002-1e-8)
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)

    def test_elastic_proposal_cannot_accept_mutually_incompatible_recipient_clearances(self):
        v=np.array([[.015,.004,.004],[.015,.006,.004],[.015,.005,.006]])
        def sheet(x,reverse=False):
            faces=np.array([[0,1,2],[0,2,3]])
            if reverse:faces=faces[:,::-1]
            return SurfaceIndex(Mesh(np.array([[x,-.01,-.01],[x,.02,-.01],[x,.02,.02],[x,-.01,.02]]),faces,'atlas-world'))
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
            SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
            recipient_surfaces=[(sheet(.0149,True),False,'oriented-exterior'),
                                (sheet(.0151),False,'oriented-exterior')],
            rigid_translation=True,elastic_recipient_proposal=True,maximum_iterations=4)
        self.assertFalse(proof['complete']);self.assertTrue(proof['finalClearanceChecksUnrelaxed'])
        self.assertTrue(all(r['feasible']for r in proof['runs']))
        self.assertTrue(any(not r['originalHalfspacesSatisfied']for r in proof['runs']))
        self.assertLess(proof['minimumSampledRecipientGap_m'],0.)
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)

    def test_soft_tissue_margin_does_not_change_bone_or_skin_margin(self):
        v=np.array([[.0144,.004,.004],[.0144,.006,.004],[.0144,.005,.006]])
        sheet=SurfaceIndex(Mesh(np.array([[.0145,-.01,-.01],[.0145,.02,-.01],[.0145,.02,.02],[.0145,-.01,.02]]),
                                np.array([[0,1,2],[0,2,3]]),'atlas-world'))
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(sheet,False,'unsigned-sheet')],rigid_translation=True,recipient_clearance_m=.00002)
        self.assertTrue(proof['complete'],proof);np.testing.assert_allclose(q,v,atol=1e-12)
        self.assertGreaterEqual(proof['minimumSampledBoneGap_m'],.0002)
        self.assertGreaterEqual(proof['minimumSampledSkinGap_m'],.0002)
        self.assertEqual(proof['recipientClearance_m'],.00002)

    def test_directed_escape_does_not_move_an_uninvolved_clear_section(self):
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006]])
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(SurfaceIndex(cube([.0135,-.01,.0045],[.0145,.02,.0055])),True,'directed-exterior','fixed rod',[1.,0.,0.])],
                    rigid_translation=True)
        self.assertTrue(proof['complete'],proof);np.testing.assert_allclose(q,v,atol=1e-12)

    def test_explicit_escape_direction_clears_a_closed_rod_with_actual_exit_rays(self):
        v=np.array([[.0133,.004,.005],[.0147,.006,.005],[.014,.005,.0058]])
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(SurfaceIndex(cube([.0135,-.01,.0045],[.0145,.02,.0055])),True,'directed-exterior','fixed rod',[0.,0.,1.])],
                    rigid_translation=True)
        self.assertTrue(proof['complete'],proof);self.assertGreater(q[:,2].min(),.00569)
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)

    def test_source_edge_patches_preserve_profile_shape_and_have_bounded_diameter(self):
        profile=np.array([[0.,0.,0.],[.002,0.,0.],[0.,.002,0.]])
        core=np.c_[np.full(6,.001),np.zeros(6),.004+np.arange(6)*.0004]
        v=np.r_[profile,core];f=np.array([[0,k,k+1]for k in range(3,8)])
        W,centres,_,_=source_control_basis(v,f,[np.arange(3)],np.arange(9),core_patch_diameter_m=.001)
        self.assertLess(len(centres),7)
        owners=W.argmax(1)
        self.assertTrue(np.all(owners[:3]==owners[0]));self.assertFalse(np.isin(owners[0],owners[3:]).any())
        for label in np.unique(owners[3:]):
            p=v[owners==label]
            self.assertLessEqual(np.linalg.norm(p[:,None]-p[None,:],axis=2).max(),.001)
        q=v+W@(np.arange(len(centres)*3).reshape(-1,3)*.0001)
        np.testing.assert_allclose(q[:3,None]-q[None,:3],v[:3,None]-v[None,:3],atol=1e-18)

    def test_nearby_unconnected_source_controls_are_not_patch_merged(self):
        v=np.array([[0.,0.,0.],[.001,0.,0.],[0.,.001,0.],[.00001,.00001,.00001]])
        W,_,_,_=source_control_basis(v,np.array([[0,1,2]]),[],np.arange(4),core_patch_diameter_m=.003)
        self.assertNotEqual(W.argmax(1)[0],W.argmax(1)[3])

    def test_one_rigid_section_uses_one_escape_plane_around_a_closed_rod(self):
        v=np.array([[.0133,.004,.005],[.0147,.006,.005],[.014,.005,.0058]])
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(SurfaceIndex(cube([.0135,-.01,.0045],[.0145,.02,.0055])),True,'coherent-section-exterior')],
                    rigid_translation=True)
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)
        self.assertTrue(np.any(q[:,0].max()<.0133 or q[:,0].min()>.0147 or q[:,2].max()<.0043 or q[:,2].min()>.0057))

    def test_harmonic_core_transport_keeps_the_source_mesh_and_reduces_variables(self):
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006],
                    [.014,.004,.004],[.014,.006,.004],[.014,.005,.006],[.013,.005,.005]])
        f=np.array([[0,1,6],[1,4,6],[4,3,6],[3,0,6],[1,2,5],[1,5,4],[0,3,5],[0,5,2]])
        profiles=[np.arange(3),np.arange(3,6)]
        sheet=SurfaceIndex(Mesh(np.array([[.0145,-.01,-.01],[.0145,.02,-.01],[.0145,.02,.02],[.0145,-.01,.02]]),
                                np.array([[0,1,2],[0,2,3]]),'atlas-world'))
        q,proof=clear_source_controls(v,f,profiles,np.arange(7),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(sheet,False,'oriented-exterior')],harmonic_cores=True,preserve_dihedral_sign=False)
        self.assertTrue(proof['complete'],proof);self.assertEqual(proof['controls'],2)
        self.assertEqual(q.shape,v.shape);self.assertLess(proof['sectionGaugeError_m'],1e-16)
        self.assertGreater(q[:,0].min(),.01469)

    def test_unsigned_open_sheet_has_no_invented_solid_interior(self):
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006]])
        sheet=SurfaceIndex(Mesh(np.array([[.0145,-.01,-.01],[.0145,.02,-.01],[.0145,.02,.02],[.0145,-.01,.02]]),
                                np.array([[0,1,2],[0,2,3]]),'atlas-world'))
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(sheet,False,'unsigned-sheet')],rigid_translation=True)
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_allclose(q,v,atol=1e-12)
        self.assertTrue(proof['actualTriangleAuditRequired'])

    def test_open_sheet_boundary_does_not_push_a_clear_source_across_its_extended_plane(self):
        v=np.array([[.012,.006,.004],[.014,.007,.004],[.014,.006,.006]])
        sheet=SurfaceIndex(Mesh(np.array([[.0145,-.005,-.01],[.0145,.005,-.01],[.0145,.005,.02],[.0145,-.005,.02]]),
                                np.array([[0,1,2],[0,2,3]]),'atlas-world'))
        q,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(sheet,False,'oriented-exterior')],rigid_translation=True)
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_allclose(q,v,atol=1e-12)

    def test_oriented_recipient_exterior_does_not_retain_an_interior_starting_side(self):
        v=np.array([[.012,.004,.004],[.014,.006,.004],[.014,.005,.006]])
        f=np.array([[0,1,2]])
        sheet=SurfaceIndex(Mesh(np.array([[.0145,-.01,-.01],[.0145,.02,-.01],[.0145,.02,.02],[.0145,-.01,.02]]),
                                np.array([[0,1,2],[0,2,3]]),'atlas-world'))
        q,proof=clear_source_controls(v,f,[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(sheet,False,'oriented-exterior')],rigid_translation=True,preserve_dihedral_sign=False)
        self.assertTrue(proof['complete'],proof);self.assertGreater(q[:,0].min(),.01469)
        self.assertFalse(proof['sourceDihedralSignPreserved'])
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)

    def test_complete_source_route_can_warm_start_without_changing_reference_calibre(self):
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006],
                    [.014,.004,.004],[.014,.006,.004],[.014,.005,.006]])
        f=np.array([[0,1,2],[3,4,5]])
        q,proof=clear_source_controls(v,f,[np.arange(3),np.arange(3,6)],np.arange(6),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),
                    SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    initial_positions=v+np.array([.001,0.,0.]))
        self.assertTrue(proof['complete'],proof)
        self.assertTrue(proof['warmStartUsed'])
        self.assertLess(proof['sectionGaugeError_m'],1e-16)

    def test_deep_closed_recipient_penetration_cannot_be_skipped_by_nearness_filter(self):
        v=np.array([[.012,.004,.004],[.014,.006,.004],[.014,.005,.006]])
        _,proof=clear_source_controls(v,np.array([[0,1,2]]),[np.arange(3)],np.arange(3),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),
                    SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(SurfaceIndex(cube([-.02,-.02,-.02],[.04,.04,.04])),True)],
                    rigid_translation=True)
        self.assertFalse(proof['complete'],proof)
        self.assertLess(proof['minimumSampledRecipientGap_m'],-.003)

    def test_capsule_policy_uses_bone_exterior_even_when_source_starts_inside_sheet(self):
        v=np.array([[.012,.004,.004],[.014,.006,.004],[.014,.005,.006],
                    [.012,.004,.007],[.014,.006,.007],[.014,.005,.009]])
        f=np.array([[0,1,2],[3,4,5]])
        # Deliberately reverse sheet winding: the fixed bone determines its
        # declared exterior; a source-side choice would retain the wrong side.
        sheet=SurfaceIndex(Mesh(np.array([[.0145,-.01,-.01],[.0145,.02,-.01],[.0145,.02,.02],[.0145,-.01,.02]]),
                                np.array([[2,1,0],[3,2,0]]),'atlas-world'))
        q,proof=clear_source_controls(v,f,[np.arange(3),np.arange(3,6)],np.arange(6),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),
                    SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(sheet,False,'bone-exterior')],rigid_translation=True)
        self.assertTrue(proof['complete'],proof)
        self.assertGreater(q[:,0].min(),.01469)
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)

    def test_rigid_translation_clears_finite_open_recipient_sheet_without_changing_source_shape(self):
        v=np.array([[.012,.004,.004],[.014,.006,.004],[.014,.005,.006],
                    [.012,.004,.007],[.014,.006,.007],[.014,.005,.009]])
        f=np.array([[0,1,2],[3,4,5]])
        sheet=SurfaceIndex(Mesh(np.array([[.0125,-.01,-.01],[.0125,.02,-.01],[.0125,.02,.02],[.0125,-.01,.02]]),
                                np.array([[0,1,2],[0,2,3]]),'atlas-world'))
        q,proof=clear_source_controls(v,f,[np.arange(3),np.arange(3,6)],np.arange(6),
                    SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])),
                    SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03])),
                    recipient_surfaces=[(sheet,False)],rigid_translation=True)
        self.assertTrue(proof['complete'],proof)
        self.assertGreater(q[:,0].min(),.01269)
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)
        self.assertEqual(proof['recipientSheetSidePolicies'],['initial-control-side'])

    def test_protected_source_mouth_stays_literal_in_corridor_solve(self):
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006],
                    [.014,.004,.004],[.014,.006,.004],[.014,.005,.006]])
        f=np.array([[0,1,2],[3,4,5]])
        bone=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01]))
        skin=SurfaceIndex(cube([-.01,-.01,-.01],[.03,.03,.03]))
        q,proof=clear_source_controls(v,f,[np.arange(3),np.arange(3,6)],np.arange(6),
                    bone,skin,protected_raw=[0,1,2])
        self.assertTrue(proof['complete'],proof)
        np.testing.assert_array_equal(q[:3],v[:3])
        self.assertEqual(proof['fixedSourceControls'],1)

    def test_section_and_aliases_keep_exact_calibre_under_any_control_translation(self):
        v=np.array([[0.,0.,0.],[.001,0.,0.],[0.,.001,0.],[0.,0.,.002],[0.,0.,0.]])
        f=np.array([[0,1,3],[1,2,3],[2,4,3]])
        W,c,_,_=source_control_basis(v,f,[np.array([0,1,2])],[0,1,2,3,0])
        delta=np.arange(len(c)*3,dtype=float).reshape(-1,3)*.0001
        q=v+W@delta
        np.testing.assert_allclose(q[:3,None]-q[None,:3],v[:3,None]-v[None,:3],atol=1e-18)
        np.testing.assert_array_equal(q[0],q[4])

    def test_micro_coherence_preserves_distinct_source_offsets_without_welding(self):
        v=np.array([[0.,0.,0.],[.000001,0.,0.],[0.,.001,0.],[0.,0.,.001]])
        f=np.array([[0,2,3],[1,3,2]])
        W,c,_,_=source_control_basis(v,f,[],np.arange(4),micro_coherence_m=.000002)
        q=v+W@(np.arange(len(c)*3).reshape(-1,3)*.0001)
        np.testing.assert_allclose(q[1]-q[0],v[1]-v[0],atol=1e-18)
        self.assertFalse(np.array_equal(q[0],q[1]))
        self.assertEqual(len(q),len(v))


if __name__=='__main__':
    unittest.main()
