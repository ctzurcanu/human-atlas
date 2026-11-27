import unittest
from copy import deepcopy
import numpy as np
from geometry_tools.algorithms.longitudinal_organ_cage import LongitudinalOrganCage,polygon_moments
from geometry_tools.algorithms.tube_sections import rotation_minimizing_frames,path_tangents
from geometry_tools.algorithms.annular_pipe_envelopes import projected_area
from geometry_tools.fit import cage_mesh
from geometry_tools.model import Mesh

class LongitudinalCageTests(unittest.TestCase):
    def test_polygon_density_and_winding_do_not_change_area_or_moments(self):
        a=np.array([[-2.,-1.],[2.,-1.],[2.,1.],[-2.,1.]])+[40.,-70.]
        b=np.insert(a,1,(a[0]+a[1])/2,axis=0)[::-1]
        ar,c,v=polygon_moments([a]);br,d,w=polygon_moments([b]);self.assertAlmostEqual(ar,8);np.testing.assert_allclose(c,[40.,-70.]);np.testing.assert_allclose(v,np.diag([4/3,1/3]));self.assertAlmostEqual(ar,br);np.testing.assert_allclose(c,d);np.testing.assert_allclose(v,w)

    def test_cross_section_orientation_does_not_flip_at_ninety_degrees(self):
        x=np.linspace(0,1,19);theta=x*np.pi*.9;source=np.broadcast_to(np.diag([4.,1.]),(len(x),2,2));r=np.array([[[np.cos(t),-np.sin(t)],[np.sin(t),np.cos(t)]]for t in theta]);target=r@source@r.transpose(0,2,1);centers=np.zeros((len(x),2));cage=LongitudinalOrganCage([0,1],[0,2],x,centers,centers,source,target)
        point=np.c_[x,np.ones(len(x)),np.zeros(len(x))];mapped=cage.map(point);np.testing.assert_allclose(mapped[:,1:],r[:,:,0],atol=1e-12);self.assertTrue(np.all(np.linalg.det(cage.jacobians(point))>0))

    def test_tube_calibre_is_independent_of_organ_flattening(self):
        x=np.linspace(0,1,25);source=np.broadcast_to(np.diag([.01,.01]),(len(x),2,2));target=np.broadcast_to(np.diag([.04,.0001]),source.shape);centers=np.c_[x*x*.02,np.sin(x)*.01];cage=LongitudinalOrganCage([0,1],[0,.85],x,np.zeros_like(centers),centers,source,target)
        axis=np.c_[x,np.zeros((len(x),2))];frames=rotation_minimizing_frames(axis,path_tangents(axis));theta=np.arange(12)*2*np.pi/12;ring=np.c_[np.cos(theta),np.sin(theta),np.zeros(12)];field=np.array([axis[:,None]+np.einsum('sij,vj->svi',frames,ring)*r for r in [.001,.0012]]);out=cage.tube_field(field,.85)
        for k in range(2):
            for i in range(len(x)):self.assertAlmostEqual(projected_area(out[k,i]),projected_area(field[k,i])*.85**2,places=16)


class ProtectedBasalLiftTests(unittest.TestCase):
    @staticmethod
    def spec():
        centers=[[0.,1.620],[.001,1.622],[0.,1.621]]
        covariances=[[[.000225,0.],[0.,.0001]]]*3
        return dict(source_bounds=[.02,.09],target_bounds=[.02,.09],fractions=[0.,.5,1.],
                    source_centers=deepcopy(centers),target_centers=deepcopy(centers),
                    source_covariances=deepcopy(covariances),target_covariances=deepcopy(covariances),
                    protected_basal_lift=dict(amplitude_m=.006,chart_axes=['Z','X','Y'],
                                             z_fall_bounds_m=[.040,.055],y_fall_bounds_m=[1.610,1.633]))

    def test_analytic_field_and_jacobian_match_independent_differences(self):
        cage=LongitudinalOrganCage(**self.spec())
        midpoint=np.array([[.0475,.002,1.6215]])
        np.testing.assert_allclose(cage.map(midpoint),midpoint+[0.,0.,.0015],atol=1e-15,rtol=0)
        expected=np.eye(3);expected[2,0]=-.3;expected[2,2]=1-.0045/.023
        np.testing.assert_allclose(cage.jacobians(midpoint)[0],expected,atol=1e-13,rtol=0)
        points=np.array([[.043,.004,1.615],[.052,-.002,1.628],[.030,0.,1.620],[.060,0.,1.615],[.0475,0.,1.640]])
        for step in [1e-7,3e-8]:
            numerical=np.stack([(cage.map(points+np.eye(3)[axis]*step)-cage.map(points-np.eye(3)[axis]*step))/(2*step)for axis in range(3)],axis=2)
            np.testing.assert_allclose(cage.jacobians(points,step=.01),numerical,atol=8e-9,rtol=0)

    def test_dorsal_and_front_protection_is_literal_and_c1(self):
        cage=LongitudinalOrganCage(**self.spec())
        protected=np.array([[.055,.003,1.610],[.2,-.001,1.600],[.04,.001,1.633],[.02,-.001,1.7]])
        np.testing.assert_array_equal(cage.map(protected),protected)
        np.testing.assert_array_equal(cage.jacobians(protected),np.broadcast_to(np.eye(3),(4,3,3)))
        # Both clipped ends of each fall have zero derivative, including the
        # fully lifted basal/posterior region, rather than a slope jump there.
        corners=np.array([[z,0.,y]for z in [.04,.055]for y in [1.610,1.633]])
        np.testing.assert_array_equal(cage.jacobians(corners),np.broadcast_to(np.eye(3),(4,3,3)))
        self.assertAlmostEqual(cage.map([[.03,0.,1.600]])[0,2],1.606)
        eps=1e-8
        near=np.array([[.04+eps,0.,1.6215],[.055-eps,0.,1.6215],[.0475,0.,1.610+eps],[.0475,0.,1.633-eps]])
        edges=np.array([[.04,0.,1.6215],[.055,0.,1.6215],[.0475,0.,1.610],[.0475,0.,1.633]])
        self.assertLess(np.max(np.abs(cage.jacobians(near)-cage.jacobians(edges))),1e-5)

    def test_global_monotone_y_bound_is_attained_and_no_other_axis_moves(self):
        cage=LongitudinalOrganCage(**self.spec())
        y=np.r_[np.linspace(1.55,1.70,1001),1.6215];y.sort()
        points=np.c_[np.full(len(y),.03),np.full(len(y),.007),y]
        mapped=cage.map(points);bound=1-1.5*.006/.023
        np.testing.assert_array_equal(mapped[:,:2],points[:,:2])
        self.assertTrue(np.all(np.diff(mapped[:,2])>0))
        det=np.linalg.det(cage.jacobians(points))
        self.assertAlmostEqual(det.min(),bound,places=13)
        self.assertAlmostEqual(cage.protected_basal_lift_jacobian_lower_bound,bound,places=13)
        self.assertTrue(np.all(det>=bound-1e-13))

    def test_rejects_nonidentity_base_and_unsupported_or_noninvertible_fields(self):
        for key,index in [('target_bounds',(1,)),('target_centers',(1,0)),('target_covariances',(1,0,0))]:
            spec=self.spec();value=np.array(spec[key]);value[index]+=1e-12;spec[key]=value
            with self.subTest(key=key),self.assertRaisesRegex(ValueError,'identical'):
                LongitudinalOrganCage(**spec)
        invalid=[('amplitude_m',-.001),('amplitude_m',.016),('amplitude_m',float('nan')),
                 ('chart_axes',['X','Y','Z']),('chart_axes',None),('chart_axes','ZXY'),
                 ('z_fall_bounds_m',[.055,.04]),('y_fall_bounds_m',[1.61,1.61])]
        for key,value in invalid:
            spec=self.spec();spec['protected_basal_lift'][key]=value
            with self.subTest(key=key,value=value),self.assertRaises(ValueError):LongitudinalOrganCage(**spec)
        spec=self.spec();spec['protected_basal_lift'].update(amplitude_m=1.,y_fall_bounds_m=[0.,1.5])
        with self.assertRaisesRegex(ValueError,'nonpositive global Jacobian'):LongitudinalOrganCage(**spec)
        spec=self.spec();spec['minimum_stretch']=1.01
        with self.assertRaisesRegex(ValueError,'identity base'):LongitudinalOrganCage(**spec)
        cage=LongitudinalOrganCage(**self.spec());cage.target_centers[0,0]+=.001
        with self.assertRaisesRegex(ValueError,'identical'):cage.map([[.04,0.,1.62]])

    def test_zero_amplitude_is_no_effect_and_cage_mesh_transports_normals(self):
        zero=self.spec();zero['protected_basal_lift']['amplitude_m']=0.
        points=np.array([[.0475,0.,1.6215],[.0485,0.,1.6215],[.0475,.001,1.6215]])
        cage=LongitudinalOrganCage(**zero)
        np.testing.assert_array_equal(cage.map(points),points)
        np.testing.assert_array_equal(cage.jacobians(points),np.broadcast_to(np.eye(3),(3,3,3)))
        attrs=dict(normals=np.tile([0,0,32767],(3,1)).astype('<i2'),
                   uvs=np.array([[0.,0.],[1.,0.],[0.,1.]],'<f4'),normalSigns=np.array([1,-1,1],'i1'))
        mesh=Mesh(points,np.array([[0,1,2]]),'rigid-Z-X-Y-test',attrs)
        result,metrics=cage_mesh(mesh,self.spec())
        # The original normal is normal to the two source tangents. Its new
        # inverse-transpose image must be perpendicular to the mapped tangents.
        normal=result.attributes['normals'][0].astype(float)/32767
        np.testing.assert_allclose(normal,np.array([.3,0.,1.])/np.sqrt(1.09),atol=1/32767,rtol=0)
        self.assertLess(abs(normal@np.array([1.,0.,-.3])),2/32767)
        self.assertEqual(normal@np.array([0.,1.,0.]),0.)
        np.testing.assert_array_equal(result.faces,mesh.faces)
        for name in ['uvs','normalSigns']:np.testing.assert_array_equal(result.attributes[name],mesh.attributes[name])
        self.assertAlmostEqual(metrics['protected_basal_lift_jacobian_lower_bound'],1-1.5*.006/.023,places=13)
        self.assertFalse(metrics['meshed_self_intersection_proof'])

if __name__=='__main__':unittest.main()
