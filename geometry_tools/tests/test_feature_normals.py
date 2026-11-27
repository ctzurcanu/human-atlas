import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.feature_normals import FeatureNormals,metric_oriented_rows,classified_oriented_distance,inward_sheet_clearance


class FeatureNormalTests(unittest.TestCase):
    def test_convex_tetrahedron_vertex_sign_does_not_depend_on_arbitrary_incident_face(self):
        v=np.array([[0.,0.,0.],[.01,0.,0.],[0.,.01,0.],[0.,0.,.01]])
        f=np.array([[0,2,1],[0,1,3],[0,3,2],[1,2,3]])
        field=FeatureNormals(Mesh(v,f,'test:metres'))
        hit=dict(face_ids=np.array([0,1,3]),barycentric=np.array([[0.,0.,1.],[0.,1.,0.],[1.,0.,0.]]),
                 points=np.repeat(v[1,None],3,axis=0))
        q=hit['points']+[.001,.0008,.0008]
        self.assertLess((q[0]-v[1])@field.face_normal[0],0.)
        normal,boundary=field.at(hit);self.assertFalse(boundary.any())
        np.testing.assert_allclose(normal,np.repeat(normal[:1],3,axis=0),atol=1e-15)
        _,gap=metric_oriented_rows(q,hit,normal);self.assertTrue((gap>0).all())
        np.testing.assert_allclose(gap,np.linalg.norm(q-hit['points'],axis=1),atol=1e-17)

    def test_open_rim_is_reported_and_does_not_extend_as_an_infinite_plane(self):
        v=np.array([[0.,0.,0.],[0.,.01,0.],[0.,.01,.01],[0.,0.,.01]])
        field=FeatureNormals(Mesh(v,np.array([[0,1,2],[0,2,3]]),'test:metres'))
        hit=dict(face_ids=np.array([0,0]),barycentric=np.array([[0.,.5,.5],[.5,0.,.5]]),points=np.array([[0.,.01,.005],[0.,.005,.005]]))
        normal,boundary=field.at(hit);np.testing.assert_array_equal(boundary,[True,False])
        q=hit['points']+[-.001,.001,0.]
        _,gap=metric_oriented_rows(q,hit,normal,unsigned_boundary=boundary)
        self.assertGreater(gap[0],0.);self.assertLess(gap[1],0.)
        signed,known=classified_oriented_distance(q,hit,normal,boundary)
        self.assertTrue(np.isnan(signed[0]));self.assertFalse(known[0])
        self.assertLess(signed[1],0.);self.assertTrue(known[1])
        self.assertEqual(int((signed>0).sum()),0)
        inward,clearance,known=inward_sheet_clearance(q,hit,normal,boundary)
        self.assertTrue((clearance>0).all())
        self.assertGreater(inward[0]@(q[0]-hit['points'][0]),0.)
        np.testing.assert_array_equal(known,[False,True])

    def test_masked_ambiguous_edge_remains_unknown_without_aborting_valid_samples(self):
        v=np.array([[0.,0.,0.],[.01,0.,0.],[0.,.01,0.]])
        field=FeatureNormals(Mesh(v,np.array([[0,1,2],[1,0,2]]),'test:metres'))
        hit=dict(face_ids=np.array([0,0]),barycentric=np.array([[.5,.5,0.],[1/3,1/3,1/3]]),points=np.array([[.005,0.,0.],[.01/3,.01/3,0.]]))
        with self.assertRaisesRegex(ValueError,'ambiguous'):field.at(hit)
        normal,unknown,valid=field.at_masked(hit)
        np.testing.assert_array_equal(valid,[False,True]);np.testing.assert_array_equal(unknown,[True,False]);np.testing.assert_array_equal(normal[0],[0.,0.,0.]);np.testing.assert_allclose(normal[1],[0.,0.,1.])
        signed,known=classified_oriented_distance(hit['points']+[0.,0.,.001],hit,normal,unknown)
        self.assertTrue(np.isnan(signed[0]));np.testing.assert_array_equal(known,[False,True]);self.assertAlmostEqual(signed[1],.001)
        strict,_=field.at({k:value[1:]for k,value in hit.items()});np.testing.assert_array_equal(strict,normal[1:])
