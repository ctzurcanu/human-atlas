import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.radial_lining import RadialLining,radial_apposition,radial_jacobians

class RadialLiningTests(unittest.TestCase):
    def cylinder(self):
        t=np.arange(24)*2*np.pi/24;p=np.r_[np.c_[np.cos(t),np.full(24,-1.),np.sin(t)],np.c_[np.cos(t),np.full(24,1.),np.sin(t)]];f=[]
        for i in range(24):
            j=(i+1)%24;f.extend([[i,j,24+j],[i,24+j,24+i]])
        f=np.asarray(f);q=p[f];n=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
        if np.sum(n[:,0]*q[:,:,0].mean(1)+n[:,2]*q[:,:,2].mean(1))>0:f=f[:,[0,2,1]]
        return Mesh(p,f,'fixture:local-metres')

    def test_real_facets_establish_radial_depth_including_shared_edges(self):
        mesh=self.cylinder();query=RadialLining(mesh);angles=np.array([0.,np.pi/24,.7,2.2]);r,w=query.radii(np.zeros(4),angles)
        self.assertAlmostEqual(r[0],1.);self.assertAlmostEqual(r[1],np.cos(np.pi/24))
        p=np.einsum('bi,bij->bj',w['barycentric'],mesh.triangles[w['face_ids']])
        np.testing.assert_allclose(p,np.c_[r*np.cos(angles),np.zeros(4),r*np.sin(angles)],atol=1e-14)

    def test_reversed_wetted_wall_is_rejected(self):
        m=self.cylinder()
        with self.assertRaises(ValueError):RadialLining(Mesh(m.vertices,m.faces[:,[0,2,1]],m.frame))

    def test_untapered_radial_sectors_preserve_area(self):
        p=np.array([[1.,0.,0.],[1.1,0.,0.],[1.2,0.,0.]])
        q,r=radial_apposition(p,np.ones(3),np.ones(3),.75,10.)
        self.assertAlmostEqual(q[0,0],.5)
        np.testing.assert_allclose(np.diff(q[:,0]**2),np.diff(p[:,0]**2),atol=1e-14)
        self.assertTrue(np.all(q[:,0]>0));self.assertTrue(np.all(np.diff(q[:,0])>0))
        self.assertTrue(r['continuous_sector_area_preserved_where_taper_is_one'])

    def test_outer_support_and_invalid_tissue_controls(self):
        p=np.array([[1.,0.,0.],[2.,0.,0.]])
        q,r=radial_apposition(p,np.ones(2),np.ones(2),.5,.2)
        np.testing.assert_array_equal(q[1],p[1])
        with self.assertRaises(ValueError):radial_apposition(np.array([[.9,0.,0.]]),np.ones(1),np.ones(1),.5,.2)


class RadialDerivativeTests(unittest.TestCase):
    def field(self,p):
        theta=np.arctan2(p[:,2],p[:,0]);r0=1+.05*p[:,1]+.02*np.sin(theta);window=.7-.1*p[:,1]
        return radial_apposition(p,r0,window,.6,1.)[0]

    def test_jacobian_matches_actual_field_and_has_positive_sector_volume(self):
        p=np.array([[1.2,.1,.25],[1.7,.2,.3]]);theta=np.arctan2(p[:,2],p[:,0]);r0=1+.05*p[:,1]+.02*np.sin(theta)
        gradients=np.c_[np.full(2,.05),.02*np.cos(theta)];j,report=radial_jacobians(p,r0,gradients,.7-.1*p[:,1],np.full(2,-.1),.6,1.)
        numeric=np.empty_like(j);step=1e-6
        for k in range(3):
            delta=np.zeros_like(p);delta[:,k]=step;numeric[:,:,k]=(self.field(p+delta)-self.field(p-delta))/(2*step)
        np.testing.assert_allclose(j,numeric,atol=1e-8,rtol=1e-8)
        np.testing.assert_allclose(np.linalg.det(j),report['continuous_determinant'],atol=1e-14)
        self.assertTrue(np.all(report['continuous_determinant']>=1.))

    def test_fully_closed_surface_cannot_reuse_a_singular_normal_transport(self):
        with self.assertRaisesRegex(ValueError,'remeshed contact'):
            radial_jacobians(np.array([[1.,0.,0.]]),np.ones(1),np.zeros((1,2)),np.ones(1),np.zeros(1),1.,1.)
