import unittest
import numpy as np
from geometry_tools.algorithms.compact_axis_field import compact_axis_field


class CompactAxisFieldTest(unittest.TestCase):
 def test_field_is_literal_outside_and_at_support_boundary(self):
  p=np.array([[0,0,0],[1,0,0],[1.1,0,0],[0,2,0]],float);r=compact_axis_field(p,[0,0,0],[1,1,1],.1)
  np.testing.assert_array_equal(r['points'][1:],p[1:]);np.testing.assert_array_equal(r['jacobians'][1:],np.repeat(np.eye(3)[None],3,axis=0));self.assertEqual(r['points'][0,2],.1)
 def test_analytic_derivative_and_bound_against_independent_difference(self):
  p=np.random.default_rng(11).uniform(-.5,.5,(15,3));radii=np.array([.8,1.2,.9]);r=compact_axis_field(p,[0,0,0],radii,.12)
  for axis in range(3):
   d=np.eye(3)[axis]*1e-6;a=compact_axis_field(p+d,[0,0,0],radii,.12)['points'];b=compact_axis_field(p-d,[0,0,0],radii,.12)['points'];np.testing.assert_allclose((a-b)/2e-6,r['jacobians'][:,:,axis],atol=1e-10)
  self.assertGreaterEqual(np.linalg.det(r['jacobians']).min(),r['continuous_determinant_lower_bound'])
 def test_inverse_transpose_normal_matches_transformed_tangent_plane(self):
  r=compact_axis_field([[.2,.1,-.3]],[0,0,0],[.8,1.2,.9],.12);J=r['jacobians'][0];a=np.array([1.,2.,0.]);b=np.array([0.,1.,2.]);n=np.cross(a,b);transported=np.linalg.solve(J.T,n);self.assertAlmostEqual(float(transported@(J@a)),0.,places=12);self.assertAlmostEqual(float(transported@(J@b)),0.,places=12)
 def test_folding_parameter_rejected(self):
  with self.assertRaisesRegex(ValueError,'positive continuous'):compact_axis_field([[0,0,0]],[0,0,0],[1,1,1],1.)
 def test_gate_keeps_the_protected_halfspace_literal_and_mirrors(self):
  p=np.array([[.1,.0,.2],[.25,.0,.2],[.5,.0,.2]]);gate=dict(axis=0,sign=1,start_m=.1,width_m=.3);a=compact_axis_field(p,[0,0,0],[1,1,1],.1,gate=gate);self.assertEqual(a['points'][0,2],p[0,2]);self.assertGreater(a['points'][1,2],p[1,2]);mirror=p.copy();mirror[:,0]*=-1;b=compact_axis_field(mirror,[0,0,0],[1,1,1],.1,gate=dict(gate,sign=-1));np.testing.assert_array_equal(a['points'][:,2],b['points'][:,2])
 def test_gated_derivative_matches_finite_difference(self):
  p=np.array([[.2,.1,-.3]]);gate=dict(axis=0,sign=1,start_m=.1,width_m=.3);r=compact_axis_field(p,[0,0,0],[.8,1.2,.9],.12,gate=gate)
  for axis in range(3):
   d=np.eye(3)[axis]*1e-6;a=compact_axis_field(p+d,[0,0,0],[.8,1.2,.9],.12,gate=gate)['points'];b=compact_axis_field(p-d,[0,0,0],[.8,1.2,.9],.12,gate=gate)['points'];np.testing.assert_allclose((a-b)/2e-6,r['jacobians'][:,:,axis],atol=1e-10)
  self.assertGreaterEqual(np.linalg.det(r['jacobians']).min(),r['continuous_determinant_lower_bound'])
