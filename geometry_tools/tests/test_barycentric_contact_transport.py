import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.barycentric_contact_transport import barycentric_contact_transport,endpoint_weight_rows,endpoint_replay,contact_translation_kernel


class BarycentricTransportTests(unittest.TestCase):
    def meshes(self):
        a=np.array([[0.,0,0],[1,0,0],[0,1,0]])
        b=np.array([[.25,.25,-1],[.25,.25,1],[.75,.25,0]])
        return [Mesh(v,np.array([[0,1,2]]),'fixture:metres')for v in [a,b]]

    def test_literal_contact_replays_on_both_original_triangles(self):
        ms=self.meshes();r=barycentric_contact_transport(ms,[dict(first=0,second=1,pairs=[[0,0]])])
        self.assertEqual(r['crossingPairs'],1);self.assertEqual(len(r['endpoints']),2)
        self.assertLess(endpoint_replay(ms,r['endpoints']),1e-12)

    def test_endpoint_transport_allows_nonrigid_profile_motion(self):
        ms=self.meshes();r=barycentric_contact_transport(ms,[dict(first=0,second=1,pairs=[[0,0]])]);w=[np.eye(6)[:3],np.eye(6)[3:]];A=endpoint_weight_rows(ms,w,r['endpoints']);_,_,v=np.linalg.svd(A,full_matrices=True);delta=v[-1]*.01
        self.assertGreater(float(np.ptp(delta)),.001)
        moved=[Mesh(m.vertices+w[i]@delta[:,None]*np.array([[0.,0.,1.]]),m.faces,m.frame)for i,m in enumerate(ms)]
        self.assertLess(endpoint_replay(moved,r['endpoints']),1e-12)

    def test_coplanar_overlap_is_kept_explicit_instead_of_fake_endpoints(self):
        m=self.meshes()[0];r=barycentric_contact_transport([m,m],[dict(first=0,second=1,pairs=[[0,0]])]);self.assertEqual(r['coplanarPairs'],1);self.assertFalse(r['endpoints']);self.assertFalse(r['anatomicalInterfaceEstablished'])

    def test_eliminated_contact_equations_keep_original_endpoint_tolerance(self):
        ms=self.meshes();r=barycentric_contact_transport(ms,[dict(first=0,second=1,pairs=[[0,0]])]);w=[np.eye(6)[:3],np.eye(6)[3:]];A=endpoint_weight_rows(ms,w,r['endpoints']);Z,p=contact_translation_kernel(A,maximum_motion_m=.002,maximum_endpoint_error_m=3e-11);d=Z@np.full((Z.shape[1],3),.0002);moved=[Mesh(m.vertices+w[i]@d,m.faces,m.frame)for i,m in enumerate(ms)]
        self.assertLess(endpoint_replay(moved,r['endpoints']),3e-11);self.assertLess(p['maximumEndpointErrorBound_m'],3e-11);self.assertLess(Z.shape[1],6)


if __name__=='__main__':unittest.main()
