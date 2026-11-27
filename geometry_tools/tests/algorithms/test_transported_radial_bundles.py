import copy,unittest
import numpy as np
from geometry_tools.algorithms.transported_radial_bundles import map_bundle,compose_bundles

def spec():
    return dict(stations=[dict(y_m=y,source_x_m=x,source_z_m=z,target_x_m=x+dx,target_z_m=z+dz,
        angle_radians=angle,core_m=core,width_m=width)for y,x,z,dx,dz,angle,core,width in [
        (.14,.10,-.030,.010,-.004,-.35,.0025,.007),
        (.25,.104,-.027,.014,-.006,-.48,.003,.008),
        (.36,.102,-.025,.011,-.005,-.39,.0027,.0075)]],
        y_core_interval_m=[.14,.36],y_support_interval_m=[.10,.40],steps=24)

class TransportedRadialBundlesTests(unittest.TestCase):
    def test_full_jacobian_in_radial_and_longitudinal_transitions(self):
        s=spec();p=np.array([[.106,.125,-.031],[.110,.225,-.023],[.105,.375,-.019]])
        _,j,_=map_bundle(p,s);h=1e-7;fd=np.empty_like(j)
        for k in range(3):
            d=np.zeros_like(p);d[:,k]=h
            fd[:,:,k]=(map_bundle(p+d,s)[0]-map_bundle(p-d,s)[0])/(2*h)
        np.testing.assert_allclose(j,fd,rtol=2e-5,atol=2e-6)
        self.assertGreater(np.max(np.abs(j[:,[0,2],1])),.01)

    def test_analytic_bound_against_dense_actual_determinants(self):
        s=spec();rng=np.random.default_rng(3301)
        p=np.column_stack([rng.uniform(.08,.135,700),rng.uniform(.101,.399,700),rng.uniform(-.05,-.005,700)])
        _,j,proof=map_bundle(p,s);det=np.linalg.det(j)
        self.assertGreater(proof['analytic_determinant_lower_bound'],0)
        self.assertGreaterEqual(float(det.min())+1e-12,proof['analytic_determinant_lower_bound'])
        np.testing.assert_array_equal(j[:,1,:],np.tile([0.,1.,0.],(len(p),1)))

    def test_exact_fixed_Y_and_outside_support(self):
        p=np.array([[0.,.22,0.],[.1,.10,-.03],[.1,.40,-.03],[.1,.09,-.03],[.1,.41,-.03]])
        q,j,_=map_bundle(p,spec());np.testing.assert_array_equal(q,p)
        np.testing.assert_array_equal(j,np.tile(np.eye(3),(len(p),1,1)))

    def test_complete_flat_section_is_one_proper_rigid_transform(self):
        s=spec()
        for r in s['stations']:
            r.update(source_x_m=.1,source_z_m=-.03,target_x_m=.112,target_z_m=-.035,angle_radians=-.43,core_m=.003,width_m=.008)
        angle=np.linspace(0,2*np.pi,24,endpoint=False)
        loops=[np.column_stack([.1+off+.0007*np.cos(angle),np.full(24,.22),-.03+.0007*np.sin(angle)])for off in [-.001,0,.001]]
        p=np.concatenate(loops);q,_,_=map_bundle(p,s)
        co,si=np.cos(-.43),np.sin(-.43);R=np.array([[co,-si],[si,co]])
        expected=(p[:,[0,2]]-[.1,-.03])@R.T+[.112,-.035]
        np.testing.assert_allclose(q[:,[0,2]],expected,atol=1e-15,rtol=0)
        np.testing.assert_allclose(np.linalg.norm(q[:,None]-q[None],axis=2),np.linalg.norm(p[:,None]-p[None],axis=2),atol=1e-15)
        for a,b in zip(loops,np.split(q,3)):
            area=lambda v:abs(np.sum(v[:,0]*np.roll(v[:,2],-1)-v[:,2]*np.roll(v[:,0],-1)))/2
            self.assertAlmostEqual(area(a),area(b),places=17)

    def test_composition_chain_rule(self):
        a=spec();b=copy.deepcopy(a)
        for r in b['stations']:
            r['source_x_m']+=.012;r['target_x_m']+=.014;r['angle_radians']*=.5
        p=np.array([[.103,.20,-.030],[.110,.28,-.024]])
        q,j1,_=map_bundle(p,a);r,j2,_=map_bundle(q,b);out,j,proof=compose_bundles(p,[a,b])
        np.testing.assert_array_equal(out,r);np.testing.assert_allclose(j,j2@j1,atol=1e-14)
        self.assertGreater(proof['analytic_determinant_lower_bound'],0)

    def test_nonpositive_or_invalid_controls_rejected(self):
        for edit in ['steps','width','interval','nan']:
            s=spec()
            if edit=='steps':s['steps']=1
            if edit=='width':s['stations'][0]['width_m']=0
            if edit=='interval':s['y_support_interval_m']=[.15,.35]
            if edit=='nan':s['stations'][1]['angle_radians']=float('nan')
            with self.assertRaises(ValueError):map_bundle(np.array([[.1,.22,-.03]]),s)

if __name__=='__main__':unittest.main()
