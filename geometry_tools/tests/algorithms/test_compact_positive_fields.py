import unittest
import numpy as np
from geometry_tools.algorithms.compact_positive_fields import compact_radial_map, compact_ellipsoidal_translation, separable_corridor_translation, compose_compact_radial_maps


class CompactPositiveFieldsTests(unittest.TestCase):
    def test_variable_course_keeps_calibre_and_has_full_analytic_jacobian(self):
        spec=dict(centerline_stations=[dict(y_m=-.2,x_m=0.,z_m=0.),dict(y_m=.2,x_m=0.,z_m=0.)],
            center_y_m=0.,core_axes_m=[.015,.10,.015],support_axes_m=[.10,.25,.10],
            y_core_interval_m=[-.1,.1],y_support_interval_m=[-.2,.2],
            translation_stations=[dict(y_m=-.2,x_m=0.,z_m=0.),dict(y_m=-.07,x_m=.004,z_m=.006),dict(y_m=.08,x_m=-.008,z_m=.002),dict(y_m=.2,x_m=0.,z_m=0.)])
        section=np.array([[-.005,.05,-.005],[.005,.05,-.005],[.005,.05,.005],[-.005,.05,.005]])
        q,_,proof=separable_corridor_translation(section,**spec)
        np.testing.assert_allclose(q-q[0],section-section[0],atol=1e-17)
        self.assertEqual(proof['kind'],'separable-corridor-course')
        p=np.random.default_rng(1286).uniform([-.11,-.22,-.11],[.11,.22,.11],(400,3))
        q,j,proof=separable_corridor_translation(p,**spec);finite=np.zeros_like(j)
        for axis in range(3):
            delta=np.eye(3)[axis]*1e-7
            finite[:,:,axis]=(separable_corridor_translation(p+delta,**spec)[0]-separable_corridor_translation(p-delta,**spec)[0])/(2e-7)
        np.testing.assert_allclose(j,finite,atol=5e-8)
        self.assertTrue(np.all(np.linalg.det(j)>=proof['analytic_determinant_lower_bound']-1e-12))
        np.testing.assert_array_equal(q[:,1],p[:,1])
        outside=(abs(p[:,0])>=.1)|(abs(p[:,2])>=.1)|(abs(p[:,1])>=.2)
        np.testing.assert_array_equal(q[outside],p[outside])
        with self.assertRaises(ValueError):separable_corridor_translation(p,**spec,translation_m=[0,0,.001])
        # Station endpoints can lie inside the support. The constant
        # continuation must join them without a derivative kink.
        inner={**spec,'translation_stations':[dict(y_m=-.07,x_m=.002,z_m=.001),dict(y_m=.08,x_m=-.004,z_m=.003)]}
        for y in [-.07,.08]:
            probes=np.array([[0,y-1e-8,0],[0,y+1e-8,0]])
            _,jac,_=separable_corridor_translation(probes,**inner)
            np.testing.assert_allclose(jac[0,:,1],jac[1,:,1],atol=2e-8)

    def test_corridor_preserves_complete_sections_through_the_Y_transition(self):
        spec=dict(centerline_stations=[dict(y_m=-.2,x_m=-.002,z_m=.001),dict(y_m=0.,x_m=0.,z_m=0.),dict(y_m=.2,x_m=.002,z_m=.001)],
            center_y_m=0.,core_axes_m=[.015,.055,.010],support_axes_m=[.08,.16,.040],translation_m=[0.,0.,-.007])
        section=np.array([[-.004,.075,-.003],[.004,.075,-.003],[.004,.075,.003],[-.004,.075,.003]])
        mapped,_,_=separable_corridor_translation(section,**spec)
        np.testing.assert_allclose(mapped-mapped[0],section-section[0],atol=1e-17)
        self.assertLess(mapped[0,2],section[0,2])
        p=np.random.default_rng(81).uniform([-.07,-.15,-.035],[.07,.15,.035],(512,3))
        _,jac,proof=separable_corridor_translation(p,**spec);finite=np.zeros_like(jac);eps=1e-7
        for axis in range(3):
            delta=np.eye(3)[axis]*eps
            finite[:,:,axis]=(separable_corridor_translation(p+delta,**spec)[0]-separable_corridor_translation(p-delta,**spec)[0])/(2*eps)
        np.testing.assert_allclose(jac,finite,atol=5e-8)
        self.assertTrue(np.all(np.linalg.det(jac)>=proof['analytic_determinant_lower_bound']-1e-12))
        outside=np.array([[.3,.075,.02],[.01,.5,.02]])
        np.testing.assert_array_equal(separable_corridor_translation(outside,**spec)[0],outside)
        with self.assertRaises(ValueError):
            separable_corridor_translation(p,**{**spec,'translation_m':[0.,.001,-.007]})
        asymmetric={**spec,'y_core_interval_m':[.03,.07],'y_support_interval_m':[-.02,.16]}
        _,aj,ap=separable_corridor_translation(p,**asymmetric)
        for axis in range(3):
            delta=np.eye(3)[axis]*eps
            finite[:,:,axis]=(separable_corridor_translation(p+delta,**asymmetric)[0]-separable_corridor_translation(p-delta,**asymmetric)[0])/(2*eps)
        np.testing.assert_allclose(aj,finite,atol=5e-8)
        self.assertTrue(np.all(np.linalg.det(aj)>=ap['analytic_determinant_lower_bound']-1e-12))
        fixed_tail=np.array([[0.,-.03,0.]])
        np.testing.assert_array_equal(separable_corridor_translation(fixed_tail,**asymmetric)[0],fixed_tail)

    def test_deep_ellipsoid_fixes_superficial_controls_and_composes(self):
        p=np.array([[0.,0.,0.],[.02,.03,.003],[0.,.02,-.032]])
        spec=dict(center=[0.,0.,0.],support_axes_m=[.08,.16,.024],
                  core_fraction=.2,translation_m=[0.,0.,-.007])
        q,j,proof=compact_ellipsoidal_translation(p,**spec)
        np.testing.assert_array_equal(q[0],p[0]+spec['translation_m'])
        np.testing.assert_array_equal(q[-1],p[-1])
        np.testing.assert_array_equal(j[-1],np.eye(3))
        self.assertTrue(np.all(np.linalg.det(j)>=proof['analytic_determinant_lower_bound']))
        fields=[dict(kind='ellipsoidal-translation',**spec),
                dict(kind='translation',center=[0.,0.,0.],core_m=.01,support_m=.2,
                     translation_m=[.0002,.0001,.0004])]
        _,jac,bound=compose_compact_radial_maps(p,fields)
        finite=np.zeros_like(jac);eps=1e-7
        for axis in range(3):
            delta=np.eye(3)[axis]*eps
            finite[:,:,axis]=(compose_compact_radial_maps(p+delta,fields)[0]-compose_compact_radial_maps(p-delta,fields)[0])/(2*eps)
        np.testing.assert_allclose(jac,finite,atol=3e-8)
        self.assertTrue(np.all(np.linalg.det(jac)>=bound['analytic_determinant_lower_bound']-1e-12))

    def test_ellipsoid_rejects_a_fold_or_invalid_support(self):
        for axes,d in [([.08,.16,.001],[0.,0.,-.007]),([.08,0.,.024],[0.,0.,-.001])]:
            with self.assertRaises(ValueError):
                compact_ellipsoidal_translation(np.zeros((1,3)),center=[0.,0.,0.],
                    support_axes_m=axes,core_fraction=.2,translation_m=d)

    def test_fixed_controls_outside_nonzero_center_remain_literal(self):
        p=np.array([[.020923712057992816,1.4095946351687114,.020685173183058698]])
        common=dict(center=[.001,1.61,-.028],core_m=.06,support_m=.15)
        for kind,extra in [('translation',{'translation_m':[.001,.002,.003]}),('twist',{'axis':[1.,2.,3.],'angle_radians':.2}),('dilation',{'growth':1.03})]:
            q,j,_=compact_radial_map(p,kind=kind,**common,**extra)
            np.testing.assert_array_equal(q,p)
            np.testing.assert_array_equal(j,np.eye(3)[None])

    def test_analytic_jacobians_match_finite_differences(self):
        p=np.array([[.014,.007,.003],[.03,-.004,.002],[0.,0.,0.],[.08,0.,0.]])
        for kind,extra in [('translation',{'translation_m':[.002,-.001,.003]}),
                           ('twist',{'axis':[1.,2.,3.],'angle_radians':.7}),
                           ('dilation',{'growth':1.15})]:
            spec=dict(center=[0.,0.,0.],core_m=.01,support_m=.06,kind=kind,**extra)
            out,j,proof=compact_radial_map(p,**spec)
            finite=np.zeros_like(j); eps=1e-7
            for k in range(3):
                d=np.eye(3)[k]*eps
                finite[:,:,k]=(compact_radial_map(p+d,**spec)[0]-compact_radial_map(p-d,**spec)[0])/(2*eps)
            np.testing.assert_allclose(j,finite,atol=2e-8)
            self.assertTrue(np.all(np.linalg.det(j)>=proof['analytic_determinant_lower_bound']-1e-10))
            np.testing.assert_array_equal(out[-1],p[-1])

    def test_radial_twist_is_exactly_volume_preserving_and_invertible(self):
        p=np.random.default_rng(4).normal(size=(100,3))*.03
        spec=dict(center=[.001,.002,.003],core_m=.01,support_m=.07,kind='twist',axis=[2.,-1.,3.])
        q,j,_=compact_radial_map(p,angle_radians=1.2,**spec)
        back,_,_=compact_radial_map(q,angle_radians=-1.2,**spec)
        np.testing.assert_allclose(back,p,atol=4e-17)
        np.testing.assert_allclose(np.linalg.det(j),1.,atol=2e-15)

    def test_folded_growth_or_unsupported_translation_is_rejected(self):
        with self.assertRaises(ValueError):
            compact_radial_map(np.zeros((1,3)),center=[0,0,0],core_m=.03,support_m=.031,
                               kind='dilation',growth=2.)
        with self.assertRaises(ValueError):
            compact_radial_map(np.zeros((1,3)),center=[0,0,0],core_m=.01,support_m=.02,
                               kind='translation',translation_m=[.1,0,0])

    def test_composed_fields_use_current_coordinates_and_chain_rule(self):
        p = np.array([[.019, .008, -.002], [.041, -.009, .004], [.12, 0., 0.]])
        fields = [
            dict(center=[0., 0., 0.], core_m=.01, support_m=.08,
                 kind='translation', translation_m=[.003, -.001, .002]),
            dict(center=[.005, .002, 0.], core_m=.015, support_m=.09,
                 kind='twist', axis=[1., 2., -1.], angle_radians=.25),
            dict(center=[0., 0., .001], core_m=.01, support_m=.1,
                 kind='dilation', growth=.98),
        ]
        q, jac, proof = compose_compact_radial_maps(p, fields)
        finite = np.zeros_like(jac)
        epsilon = 1e-7
        for axis in range(3):
            delta = np.eye(3)[axis] * epsilon
            finite[:, :, axis] = (
                compose_compact_radial_maps(p + delta, fields)[0]
                - compose_compact_radial_maps(p - delta, fields)[0]
            ) / (2 * epsilon)
        np.testing.assert_allclose(jac, finite, atol=3e-8)
        self.assertTrue(np.all(np.linalg.det(jac) >= proof['analytic_determinant_lower_bound'] - 1e-12))
        np.testing.assert_array_equal(q[-1], p[-1])
        sequential = p
        for field in fields:
            sequential = compact_radial_map(sequential, **field)[0]
        np.testing.assert_array_equal(q, sequential)

    def test_composition_rejects_an_uncertified_member(self):
        with self.assertRaises(ValueError):
            compose_compact_radial_maps(np.zeros((1, 3)), [
                dict(center=[0., 0., 0.], core_m=.01, support_m=.011,
                     kind='translation', translation_m=[.1, 0., 0.])
            ])

if __name__=='__main__':unittest.main()
