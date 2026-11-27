import unittest
import numpy as np
from geometry_tools.algorithms.faithful_float32 import directed_rounding
from geometry_tools.algorithms.tube_sections import TubeGeometryError
from geometry_tools.distance import triangle_pairs


class FaithfulFloat32(unittest.TestCase):
    def fixture(self):
        p=np.array([[0,0,1],[1,0,1],[0,1,1],[.1,.1,1+2e-8],[.8,.1,1+2e-8],[.1,.8,1+2e-8]],float)
        return p,p.astype(np.float32).astype(float),np.array([[0,1,2],[3,4,5]])

    def test_distinct_source_sheets_survive_encoding_without_source_motion(self):
        p,q,f=self.fixture();before=p.copy();self.assertEqual(triangle_pairs(q[f[:1]],q[f[1:]])['relation'][0],'coplanar-overlap')
        result,evidence=directed_rounding(p,q,f,f[None],time_limit=3)
        self.assertTrue(np.array_equal(p,before));self.assertEqual(triangle_pairs(result[f[:1]],result[f[1:]])['relation'][0],'near')
        self.assertTrue(evidence['faithful_floor_ceil_only']);self.assertEqual(evidence['known_crossings_remaining'],0)

    def test_protected_controls_are_not_moved_to_hide_a_failure(self):
        p,q,f=self.fixture()
        with self.assertRaises(TubeGeometryError):directed_rounding(p,q,f,f[None],mutable_controls=[],time_limit=3)

    def test_an_intersecting_anatomical_source_is_rejected(self):
        p,q,f=self.fixture();p[3:,2]=1
        with self.assertRaisesRegex(TubeGeometryError,'source already intersects'):directed_rounding(p,q,f,f[None],time_limit=3)

    def test_authored_corner_fields_survive_separation(self):
        p,q,f=self.fixture()
        words=np.broadcast_to(np.array([[1000,0,30000],[0,1000,30000],[-1000,0,30000]],np.int16),(2,3,3)).copy()
        result,_=directed_rounding(p,q,f,f[None],corner_words=words,time_limit=3)
        t=result[f];normal=np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0])
        self.assertTrue(np.all(np.einsum('fci,fi->fc',words,normal)>0))
        self.assertEqual(triangle_pairs(result[f[:1]],result[f[1:]])['relation'][0],'near')

    def test_deadline_does_not_restart_for_another_solver_round(self):
        p,q,f=self.fixture()
        with self.assertRaises(TimeoutError):directed_rounding(p,q,f,f[None],time_limit=0)

    def test_coplanar_contact_uses_a_separating_feature_in_one_round(self):
        p=np.array([[0,0,0],[1,0,0],[0,1,0],[1,1,0],[.1,.900000001,0]],float)
        f=np.array([[0,1,2],[1,3,4]])
        q=p.astype(np.float32).astype(float)
        self.assertEqual(triangle_pairs(q[f[:1]],q[f[1:]])['relation'][0],'coplanar-overlap')
        result,evidence=directed_rounding(p,q,f,f[None],mutable_controls=[4],maximum_rounds=1,time_limit=3)
        self.assertNotEqual(triangle_pairs(result[f[:1]],result[f[1:]])['relation'][0],'coplanar-overlap')
        np.testing.assert_array_equal(result[:4],q[:4])
        self.assertEqual(evidence['solver_rounds'],1)

    def test_orientation_only_recovers_a_collapsed_facet_with_fixed_native_controls(self):
        p=np.array([[0,1,0],[1,1,0],[.5,1+2e-8,0]],float)
        q=p.astype(np.float32).astype(float);f=np.array([[0,1,2]])
        words=np.tile(np.array([0,0,32767],np.int16),(1,3,1))
        self.assertEqual(np.linalg.norm(np.cross(q[1]-q[0],q[2]-q[0])),0)
        source=p.copy()
        result,evidence=directed_rounding(p,q,f,None,corner_words=words,
            mutable_controls=[2],additional_controls=[2],time_limit=3)
        np.testing.assert_array_equal(p,source)
        np.testing.assert_array_equal(result[:2],q[:2])
        self.assertGreater(np.cross(result[1]-result[0],result[2]-result[0])[2],0)
        self.assertTrue(evidence['orientation_only'])
        self.assertEqual(evidence['known_pairs_checked'],0)

    def test_orientation_only_requires_a_complete_explicit_diagnosis(self):
        p,q,f=self.fixture();words=np.tile(np.array([0,0,32767],np.int16),(2,3,1))
        with self.assertRaisesRegex(TubeGeometryError,'explicitly diagnosed'):
            directed_rounding(p,q,f,None,corner_words=words,time_limit=3)
        with self.assertRaisesRegex(TubeGeometryError,'corner fields'):
            directed_rounding(p,q,f,None,additional_controls=[3],time_limit=3)
        with self.assertRaisesRegex(TubeGeometryError,'actual indexed triangle pairs'):
            directed_rounding(p,q,f,np.empty((0,2,3),int),corner_words=words,time_limit=3)

    def test_nearby_sheet_is_retained_when_the_known_pair_is_corrected(self):
        xy=np.array([[0,0],[1,0],[0,1]],float)
        p=np.concatenate([np.c_[xy,np.full(3,z)]for z in (1-2e-8,1+2e-8,1+4e-8)])
        q=p.astype(np.float32).astype(float);f=np.arange(9).reshape(3,3)
        known=f[[0,1]][None];near=f[[1,2]][None]
        result,evidence=directed_rounding(p,q,f,known,mutable_controls=np.arange(9),
            validation_pairs=near,maximum_rounds=32,time_limit=3)
        for pair in (known,near):
            self.assertNotIn(triangle_pairs(result[pair[:,0]],result[pair[:,1]])['relation'][0],
                ('crossing','coplanar-overlap'))
        self.assertEqual(evidence['validation_pairs_checked'],1)
        self.assertTrue(evidence['source_positions_literal'])

    def test_validation_neighborhood_cannot_hide_a_source_intersection(self):
        p,q,f=self.fixture()
        with self.assertRaisesRegex(TubeGeometryError,'source already intersects'):
            directed_rounding(p,q,f,f[None],validation_pairs=np.array([[f[0],f[0]]]),time_limit=3)


if __name__=='__main__':unittest.main()
