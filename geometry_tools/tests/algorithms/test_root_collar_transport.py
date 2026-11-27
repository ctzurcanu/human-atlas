import unittest
import numpy as np
from scipy.spatial.transform import Rotation
from geometry_tools.model import Mesh
from geometry_tools.audit import audit
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.algorithms.root_collar_transport import subdivide_root_collar,transport_root_collar


def finite_pipe(n=8):
    angle=np.arange(n)*2*np.pi/n;circle=np.c_[np.cos(angle),np.sin(angle)]
    p=np.concatenate([np.c_[radius*circle,np.full(n,z)]for radius,z in [(1.,0.),(1.,1.),(.8,0.),(.8,1.)]])
    f=[]
    for i,j in zip(range(n),np.roll(np.arange(n),-1)):
        a,b,c,d=i,j,n+i,n+j
        f.extend([[a,b,d],[a,d,c]])
        a,b,c,d=2*n+i,2*n+j,3*n+i,3*n+j
        f.extend([[a,d,b],[a,c,d]])
        f.extend([[i,2*n+j,j],[i,2*n+i,2*n+j]])
        f.extend([[n+i,n+j,3*n+j],[n+i,3*n+j,3*n+i]])
    return p,np.array(f),np.r_[np.arange(n),np.arange(2*n,3*n)]


class RootCollarTransportTests(unittest.TestCase):
    def test_source_subdivision_preserves_surface_area_topology_and_open_lumen(self):
        p,f,root=finite_pipe();out=subdivide_root_collar(p,f,root)
        np.testing.assert_array_equal(out['vertices'][:len(p)],p)
        old=p[f];new=out['vertices'][out['faces']]
        area=lambda t:np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
        replay=np.bincount(out['source_face_ids'],weights=area(new),minlength=len(f))
        np.testing.assert_allclose(replay,area(old),rtol=1e-13,atol=1e-14)
        result=audit(Mesh(out['vertices'],out['faces'],'fixture:metres'),require_closed=True,expected_components=1,expected_genus=[1])
        self.assertTrue(result['passed'],result['failures'])

    def test_every_paired_profile_is_rigid_and_the_root_is_literal(self):
        p,f,root=finite_pipe();out=subdivide_root_collar(p,f,root)
        r=Rotation.from_euler('y',10,degrees=True).as_matrix();t=np.array([.05,.02,0.])
        q=transport_root_collar(out['vertices'],out['station'],r,t,[0.,0.,0.])
        np.testing.assert_array_equal(q[root],p[root])
        far=out['station']==1
        np.testing.assert_allclose(q[far],out['vertices'][far]@r.T+t,rtol=0,atol=1e-14)
        for profile in out['profiles']:
            for ring in profile['rings']:
                self.assertAlmostEqual(projected_area(q[ring]),projected_area(out['vertices'][ring]),places=13)
                before=np.linalg.norm(out['vertices'][ring,None]-out['vertices'][ring][None],axis=2)
                after=np.linalg.norm(q[ring,None]-q[ring][None],axis=2)
                np.testing.assert_allclose(after,before,rtol=0,atol=2e-14)

    def test_no_reflection_or_invalid_station_can_be_hidden(self):
        p,f,root=finite_pipe();out=subdivide_root_collar(p,f,root)
        with self.assertRaisesRegex(ValueError,'proper rigid'):
            transport_root_collar(out['vertices'],out['station'],np.diag([-1.,1.,1.]),[0.,0.,0.],[0.,0.,0.])
        with self.assertRaisesRegex(ValueError,'stations'):
            subdivide_root_collar(p,f,root,levels=[.5,.25])
        with self.assertRaisesRegex(ValueError,'Distinct'):
            subdivide_root_collar(p,f,np.r_[root,root[0]])

    def test_small_bend_remains_a_single_finite_wall(self):
        p,f,root=finite_pipe();out=subdivide_root_collar(p,f,root)
        q=transport_root_collar(out['vertices'],out['station'],Rotation.from_euler('y',8,degrees=True).as_matrix(),[.03,0.,0.],[0.,0.,0.])
        result=audit(Mesh(q,out['faces'],'fixture:metres'),require_closed=True,expected_components=1,expected_genus=[1])
        self.assertTrue(result['passed'],result['failures'])


if __name__=='__main__':unittest.main()
