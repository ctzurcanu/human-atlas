import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import SurfaceIndex
from geometry_tools.algorithms.rigid_section_clearance import fit_section_clearance


class RigidSectionClearanceTests(unittest.TestCase):
    def test_hard_skin_plane_limits_a_conflicting_clearance_objective(self):
        body=Mesh([[0,-1,-1],[0,1,-1],[0,0,1]],[[0,1,2]],'test:metres')
        v=np.array([[-.001,-.001,0],[-.001,.001,0],[-.001,0,.001],[1,1,1]])
        q,p=fit_section_clearance(v,[[0,1,2]],[np.arange(3),np.array([3])],[(SurfaceIndex(body),1,'opposing-tissue')],
            protected=[3],maximum_steps=1,wall_seconds=3,trial_scales=(1.,),
            direction_planes=dict(owner=[[0,1]],normal=[[1.,0,0]],rhs=[-.00025]))
        self.assertTrue(np.all(q[:3,0]<=-.00075+2e-9))
        self.assertTrue(np.all(q[:3,0]>-.001))
        np.testing.assert_array_equal(q[3],v[3])
        np.testing.assert_allclose(q[:3]-v[:3],np.tile(q[0]-v[0],(3,1)),atol=1e-15,rtol=0)

    def test_clearance_gain_cannot_buy_a_new_intrinsic_crossing(self):
        # Both source cells remain rigid and keep their face normals. The
        # proposed translation improves the target-plane clearance but sends
        # one cell through the other: section calibre/orientation alone pass.
        from geometry_tools.contact import self_contact_summary
        body=Mesh([[-.001,-1,-1],[-.001,1,-1],[-.001,0,1]],[[0,1,2]],'test:metres')
        v=np.array([[-.002,-.002,0],[.002,-.002,0],[0,.002,0],
                    [-.0025,-.001,-.001],[-.0025,.001,.001],[-.0025,.001,-.001]])
        faces=np.array([[0,1,2],[3,4,5]])
        saved=dict(vertices=v,owner=np.array([0,0,0,1,1,1]),shift=np.array([[0.,0,0],[.003,0,0]]))
        dangerous=v+saved['shift'][saved['owner']]
        audit,_=self_contact_summary(Mesh(dangerous,faces,'test:metres'))
        self.assertEqual(audit['defect_pairs'],1)
        q,proof=fit_section_clearance(v,faces,[np.arange(3),np.arange(3,6)],
            [(SurfaceIndex(body),1,'bone-side')],protected=[0],step_m=.003,
            maximum_steps=1,wall_seconds=3,prepared_step=saved,
            trial_scales=(1.,),reference_vertices=v,face_system=[0,0])
        np.testing.assert_array_equal(q,v)
        step=proof['steps'][0]
        self.assertEqual(step['trials'][0]['newCrossingClasses']['arterial'],1)
        self.assertLess(step['trials'][0]['energy_m2'],step['energy_m2'])
        self.assertNotIn('acceptedScale',step)

    def test_saved_step_is_bound_to_pose_and_replays_separately(self):
        from unittest.mock import patch
        from geometry_tools.contact import self_contact_summary
        body=Mesh([[0,-1,-1],[0,1,-1],[0,0,1]],[[0,1,2]],'test:metres')
        v=np.array([[-.001,-.001,0],[-.001,.001,0],[-.001,0,.001]])
        saved={};surfaces=[(SurfaceIndex(body),1,'bone-side')]
        original,p=fit_section_clearance(v,[[0,1,2]],[np.arange(3)],surfaces,
            maximum_steps=1,wall_seconds=3,prepare_step=lambda values:saved.update(values),reference_vertices=v)
        np.testing.assert_array_equal(original,v)
        with patch('geometry_tools.contact.self_contact_summary',wraps=self_contact_summary)as audit:
            q,p=fit_section_clearance(v,[[0,1,2]],[np.arange(3)],surfaces,
                maximum_steps=1,wall_seconds=3,prepared_step=saved,trial_scales=(1.,),reference_vertices=v)
            self.assertEqual(audit.call_count,1)  # The new geometry is audited; source/base caches replay.
        self.assertEqual(p['final'][0]['movableViolations'],0)
        with self.assertRaises(ValueError):
            fit_section_clearance(v+[.0001,0,0],[[0,1,2]],[np.arange(3)],surfaces,
                maximum_steps=1,wall_seconds=3,prepared_step=saved,trial_scales=(1.,),reference_vertices=v)

    def test_actual_crossing_clears_with_fixed_other_cell(self):
        from geometry_tools.contact import self_contact_summary
        v=np.array([[-.002,-.002,0],[.002,-.002,0],[0,.002,0],
                    [0,-.001,-.001],[0,.001,.001],[0,.001,-.001]])
        f=np.array([[0,1,2],[3,4,5]])
        source=v.copy();source[3:,0]+=.004
        q,p=fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],protected=[0],
            maximum_steps=4,wall_seconds=3,reference_vertices=source)
        np.testing.assert_array_equal(q[:3],v[:3])
        audit,_=self_contact_summary(Mesh(q,f,'test:metres'))
        self.assertEqual(audit['defect_pairs'],0)
        self.assertEqual(p['originalSourceContactPairsExplicitlyRetained'],0)

    def test_member_guard_counts_existing_mixed_contacts_and_keeps_fixed_member(self):
        from geometry_tools.contact import self_contact_summary
        v=np.array([[-.002,-.002,0],[.002,-.002,0],[0,.002,0],
                    [0,-.001,-.001],[0,.001,.001],[0,.001,-.001]])
        f=np.array([[0,1,2],[3,4,5]])
        q,p=fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],protected=[0],
            maximum_steps=4,wall_seconds=3,reference_vertices=v,face_system=[0,1],face_member=[0,1],resolve_mixed_contacts=True)
        self.assertEqual(p['steps'][0]['newCrossingClasses']['mixed'],1)
        self.assertEqual(p['steps'][0]['newCrossingClasses']['intrinsic'],0)
        self.assertEqual(p['steps'][0]['newCrossingClasses']['newMemberPairs'],0)
        np.testing.assert_array_equal(q[:3],v[:3]);audit,_=self_contact_summary(Mesh(q,f,'test:metres'));self.assertEqual(audit['defect_pairs'],0)

    def test_staged_member_guard_retains_all_counters_and_rejects_changed_ownership(self):
        from unittest.mock import patch
        from geometry_tools.contact import self_contact_summary
        v=np.array([[-.002,-.002,0],[.002,-.002,0],[0,.002,0],
                    [0,-.001,-.001],[0,.001,.001],[0,.001,-.001]])
        f=np.array([[0,1,2],[3,4,5]])
        saved={};kwargs=dict(protected=[0],maximum_steps=1,wall_seconds=3,
            reference_vertices=v,face_system=[0,1],face_member=[0,1],resolve_mixed_contacts=True)
        original,first=fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],
            prepare_step=lambda x:saved.update(x),**kwargs)
        np.testing.assert_array_equal(original,v)
        with patch('geometry_tools.contact.self_contact_summary',wraps=self_contact_summary)as audit:
            q,proof=fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],prepared_step=saved,**kwargs)
            self.assertGreaterEqual(audit.call_count,1)
        self.assertEqual(proof['steps'][0]['newCrossingClasses'],first['steps'][0]['newCrossingClasses'])
        self.assertEqual(proof['steps'][0]['newCrossingClasses']['mixed'],1)
        self.assertEqual(proof['steps'][0]['newCrossingClasses']['intrinsic'],0)
        self.assertEqual(proof['steps'][0]['newCrossingClasses']['newMemberPairs'],0)
        direct,direct_proof=fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],**kwargs)
        np.testing.assert_allclose(q,direct,atol=1e-14,rtol=0)
        self.assertEqual(proof['steps'][0]['trials'],direct_proof['steps'][0]['trials'])
        self.assertLessEqual(self_contact_summary(Mesh(q,f,'test:metres'))[0]['defect_pairs'],1)
        with self.assertRaises(ValueError):
            fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],prepared_step=saved,**dict(kwargs,face_member=[0,0]))

    def test_placement_only_retains_explicit_defects_without_removing_the_guards(self):
        v=np.array([[-.002,-.002,0],[.002,-.002,0],[0,.002,0],
                    [0,-.001,-.001],[0,.001,.001],[0,.001,-.001]])
        f=np.array([[0,1,2],[3,4,5]])
        saved={}
        q,proof=fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],maximum_steps=1,wall_seconds=3,
            reference_vertices=v,face_system=[0,1],face_member=[0,1],resolve_existing_member_contacts=False,
            prepare_step=lambda x:saved.update(x))
        np.testing.assert_array_equal(q,v)
        counts=proof['steps'][0]['newCrossingClasses']
        self.assertEqual(counts['mixed'],1);self.assertEqual(counts['intrinsic'],0);self.assertEqual(counts['newMemberPairs'],0)
        self.assertEqual(len(saved['base_rhs']),0)
        self.assertFalse(proof['resolveExistingMemberContacts'])
        with self.assertRaises(ValueError):
            fit_section_clearance(v,f,[np.arange(3),np.arange(3,6)],[],maximum_steps=1,wall_seconds=3,
                reference_vertices=v,face_system=[0,1],face_member=[0,1],prepared_step=saved)

    def test_actual_plane_translation_preserves_source_cell(self):
        body=Mesh([[0,-1,-1],[0,1,-1],[0,0,1]],[[0,1,2]],'test:metres')
        v=np.array([[-.001,-.001,0],[-.001,.001,0],[-.001,0,.001]])
        q,p=fit_section_clearance(v,[[0,1,2]],[np.arange(3)],[(SurfaceIndex(body),1,'bone-side')],maximum_steps=3,wall_seconds=3)
        self.assertTrue(p['constructionComplete']);self.assertEqual(p['final'][0]['movableViolations'],0)
        np.testing.assert_allclose(q[:,1:],v[:,1:],atol=1e-14)
        np.testing.assert_allclose(np.linalg.norm(q[:,None]-q[None],axis=2),np.linalg.norm(v[:,None]-v[None],axis=2),atol=1e-14)
        np.testing.assert_array_equal(body.vertices,[[0,-1,-1],[0,1,-1],[0,0,1]])

    def test_protected_cell_is_literal_and_violation_is_reported(self):
        body=Mesh([[0,-1,-1],[0,1,-1],[0,0,1]],[[0,1,2]],'test:metres')
        v=np.array([[-.001,-.001,0],[-.001,.001,0],[-.001,0,.001]])
        q,p=fit_section_clearance(v,[[0,1,2]],[np.arange(3)],[(SurfaceIndex(body),1,'bone-side')],protected=[0],maximum_steps=3,wall_seconds=3)
        np.testing.assert_array_equal(q,v)
        self.assertGreater(p['final'][0]['protectedViolations'],0)
        self.assertTrue(p['fixedControlsLiteral'])


if __name__=='__main__':unittest.main()
