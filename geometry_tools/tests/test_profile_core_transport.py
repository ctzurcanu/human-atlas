import unittest
import numpy as np
from geometry_tools.algorithms.profile_core_transport import transport_profile_cores
def tube():
 ring=np.array([[-1,-1],[1,-1],[1,1],[-1,1]],float);v=np.concatenate([np.c_[ring,np.full(4,z)]for z in range(3)]);f=[]
 for j in range(2):
  for k in range(4):a=j*4+k;b=j*4+(k+1)%4;c=a+4;d=b+4;f.extend([[a,b,d],[a,d,c]])
 return v,np.array(f)
class ProfileCoreTransportTest(unittest.TestCase):
 def test_core_interpolates_and_profiles_remain_rigid(self):
  v,f=tube();r=transport_profile_cores(v,f,[np.arange(4),np.arange(8,12)],[[0,0,0],[1,0,2]],np.repeat(np.eye(3)[None],2,axis=0));np.testing.assert_allclose(r['vertices'][4:8]-v[4:8],np.repeat([[.5,0,0]],4,axis=0),atol=1e-12);np.testing.assert_array_equal(r['vertices'][:4],v[:4]);np.testing.assert_allclose(r['vertices'][8:]-v[8:],np.repeat([[1,0,0]],4,axis=0));self.assertLess(r['harmonic_equation_maximum_residual_m'],1e-12)
 def test_explicit_source_alias_offset_rotates_without_welding(self):
  v,f=tube();v=np.r_[v,v[:1]+[1e-9,0,0]];ids=np.r_[np.arange(12),0];R=np.array([[0,-1,0],[1,0,0],[0,0,1]],float);r=transport_profile_cores(v,f,[np.arange(4),np.arange(8,12)],[[0,0,0],[0,0,2]],np.array([R,np.eye(3)]),physical_ids=ids);np.testing.assert_allclose(r['vertices'][12]-r['vertices'][0],R@(v[12]-v[0]),atol=2e-16);self.assertGreater(np.linalg.norm(r['vertices'][12]-r['vertices'][0]),0)
 def test_protected_profile_is_not_silently_released(self):
  v,f=tube()
  with self.assertRaisesRegex(ValueError,'protected'):transport_profile_cores(v,f,[np.arange(4)],[[1,0,0]],[np.eye(3)],fixed_vertices=[0])
 def test_unanchored_component_keeps_actual_source(self):
  v,f=tube();extra=np.array([[5,0,0],[6,0,0],[5,1,0]],float);v=np.r_[v,extra];f=np.r_[f,[[12,13,14]]];r=transport_profile_cores(v,f,[np.arange(4),np.arange(8,12)],[[0,0,0],[1,0,2]],np.repeat(np.eye(3)[None],2,axis=0));np.testing.assert_array_equal(r['vertices'][12:],extra);self.assertEqual(len(r['unanchored_components_retained_literal']),1)
 def test_explicit_nonprofile_port_is_retained_by_the_core_solve(self):
  v,f=tube();target=v[4]+[.25,0,0]
  r=transport_profile_cores(v,f,[np.arange(4),np.arange(8,12)],[[0,0,0],[1,0,2]],np.repeat(np.eye(3)[None],2,axis=0),prescribed_vertices=[4],prescribed_positions=[target])
  np.testing.assert_array_equal(r['vertices'][4],target);self.assertLess(r['harmonic_equation_maximum_residual_m'],1e-12)
 def test_prescribed_control_cannot_override_a_fixed_point(self):
  v,f=tube()
  with self.assertRaisesRegex(ValueError,'conflicts'):
   transport_profile_cores(v,f,[np.arange(4)],[[0,0,0]],[np.eye(3)],fixed_vertices=[4],prescribed_vertices=[4],prescribed_positions=[v[4]+1])
if __name__=='__main__':unittest.main()
