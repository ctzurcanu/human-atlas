import unittest
import numpy as np
from geometry_tools.algorithms.paired_bone_transport import paired_bone_transport,surface_frame_carriage,paired_vascular_frame_carriage,paired_curve_carriage


class PairedBoneTransportTests(unittest.TestCase):
    def test_curved_pair_maps_stations_while_retaining_source_section_calibre(self):
        t=np.linspace(0,.04,9);rails=[np.c_[np.full(9,x),np.zeros(9),t]for x in [-.004,.004]]
        bend=.006*np.sin(t/.04*np.pi);target=[r+np.c_[np.zeros(9),bend,np.zeros(9)]for r in rails]
        offsets=np.array([[0.,-.0003,0.],[0.,.00015,.00026],[0.,.00015,-.00026]])
        centres=np.array([rails[0][2],rails[1][6]])
        v=np.concatenate([c+offsets for c in centres]);sections=[np.arange(3),np.arange(3,6)]
        q,R,_=paired_curve_carriage(v,sections,np.arange(6),rails,target,[0,0,1],[0,0,1],[0,0,0],[0,0,0],1.)
        for s,expected in zip(sections,[target[0][2],target[1][6]]):
            np.testing.assert_allclose(q[s].mean(0),expected,atol=1e-16)
            np.testing.assert_allclose(np.linalg.norm(q[s,None]-q[None,s],axis=2),np.linalg.norm(v[s,None]-v[None,s],axis=2),atol=1e-16)
        np.testing.assert_allclose(np.linalg.det(R),1.,atol=1e-12)

    def test_crossing_or_coincident_rails_are_rejected(self):
        c=np.array([[0.,0.,0.],[0.,0.,.01],[0.,0.,.02]])
        with self.assertRaisesRegex(ValueError,'separated'):
            paired_curve_carriage(c,[],np.arange(3),[c,c],[c,c],[0,0,1],[0,0,1],[0,0,0],[0,0,0],1.)

    def test_both_vascular_rails_match_without_scaling_their_source_sections(self):
        source=np.array([[-.003,0.,0.],[.003,0.,0.]]);target=np.array([[.02,-.004,0.],[.02,.004,0.]])
        v=np.r_[source[0]+np.array([[0.,0.,0.],[0.,.0003,0.],[0.,0.,.0003]]),
                source[1]+np.array([[0.,0.,0.],[0.,.0003,0.],[0.,0.,.0003]])]
        sections=[np.arange(3),np.arange(3,6)];centres=np.array([v[s].mean(0)for s in sections])
        q,_,proof=paired_vascular_frame_carriage(v,sections,np.arange(6),centres,target,[0,0,1],[0,0,1],1.2)
        for s,point in zip(sections,target):
            np.testing.assert_allclose(q[s].mean(0),point,atol=1e-17)
            np.testing.assert_allclose(np.linalg.norm(q[s,None]-q[None,s],axis=2),np.linalg.norm(v[s,None]-v[None,s],axis=2),atol=1e-17)
        self.assertGreater(proof['carrierDeterminant'],0.)

    def test_surface_carrier_scales_position_and_retains_literal_section_shape(self):
        from geometry_tools.open3d_backend import SurfaceIndex
        from geometry_tools.tests.test_sparse_group_clearance import cube
        v=np.array([[.012,.004,.004],[.012,.006,.004],[.012,.005,.006]])
        q,R,proof=surface_frame_carriage(v,[np.arange(3)],[dict(sourceFrame=np.eye(3),targetFrame=np.eye(3),
                sourceCentre=np.zeros(3),targetCentre=np.array([.02,0.,0.]),scales=[1.5,1.5,1.4],
                index=SurfaceIndex(cube([0.,0.,0.],[.01,.01,.01])))],np.arange(3))
        np.testing.assert_allclose(q.mean(0),v.mean(0)*[1.5,1.5,1.4]+[.02,0.,0.],atol=1e-17)
        np.testing.assert_allclose(q-q[0],v-v[0],atol=1e-17)
        self.assertTrue(proof['sourceSectionCalibreUnscaled'])

    def test_measured_bone_width_changes_centre_offset_but_preserves_vessel_section(self):
        a=np.array([[0.,0.,0.],[0.,.02,0.],[0.,.04,0.]])
        v=np.array([[.002,.01,0.],[.003,.01,0.],[.002,.011,.001]])
        q,_,proof=paired_bone_transport(v,[np.arange(3)],a,a,[0.,0.,1.],[0.,0.,1.],
                                       transverse_station_scales=np.full((3,2),1.5))
        np.testing.assert_allclose(np.linalg.norm(q[:,None]-q[None,:],axis=2),np.linalg.norm(v[:,None]-v[None,:],axis=2),atol=1e-17)
        self.assertGreater(q.mean(0)[0],v.mean(0)[0]*1.49)
        self.assertTrue(proof['widthScalingActsOnSectionCentresOnly'])

    def test_identity_is_literal_and_source_section_calibre_survives_curved_recipient(self):
        a=np.array([[0.,0.,0.],[0.,.02,0.],[0.,.04,0.]])
        v=np.array([[.002,.01,0.],[.003,.01,0.],[.002,.011,.001],
                    [.002,.03,0.],[.003,.03,0.],[.002,.031,.001]])
        sections=[np.arange(3),np.arange(3,6)];normal=np.array([0.,0.,1.])
        q,_,_=paired_bone_transport(v,sections,a,a,normal,normal)
        np.testing.assert_array_equal(q,v)
        b=a.copy();b[1:,0]=[.006,.014]
        q,R,_=paired_bone_transport(v,sections,a,b,normal,normal)
        for s in sections:
            np.testing.assert_allclose(np.linalg.norm(q[s,None]-q[None,s],axis=2),
                                       np.linalg.norm(v[s,None]-v[None,s],axis=2),atol=1e-17)
        np.testing.assert_allclose(np.linalg.det(R),1.,atol=1e-12)

    def test_missing_or_repeated_bone_stations_are_rejected(self):
        with self.assertRaisesRegex(ValueError,'Distinct'):
            paired_bone_transport(np.zeros((3,3)),[],np.zeros((2,3)),np.zeros((2,3)),[0,0,1],[0,0,1])
