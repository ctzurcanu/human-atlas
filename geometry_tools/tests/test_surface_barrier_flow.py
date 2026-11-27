import unittest
import numpy as np
from geometry_tools.algorithms.compact_translation import CompactTranslation
from geometry_tools.algorithms.surface_barrier_flow import TriangleSurfaceGuard,guarded_translation_flow

class BarrierFlowTests(unittest.TestCase):
 def test_rigid_ellipsoidal_core_and_branch_controls_share_one_guarded_map(self):
  from geometry_tools.algorithms.surface_barrier_flow import EllipsoidalTranslationField
  field=EllipsoidalTranslationField([0,0,0],[.04,.10,.10],.8,[.00055,-.00025,-.0002])
  x=np.array([[.001,.002,.003],[.035,.02,.01],[.041,0,0],[.001,.002,.003]])
  q,J=field.map_and_jacobians(x)
  np.testing.assert_allclose(q[0]-x[0],field.displacement,atol=1e-17)
  np.testing.assert_array_equal(q[2],x[2]);np.testing.assert_array_equal(q[0],q[3])
  columns=[]
  for axis in range(3):
   e=np.eye(3)[axis]*1e-8
   columns.append((field.map(x+e)-field.map(x-e))/(2e-8))
  np.testing.assert_allclose(J,np.stack(columns,axis=2),rtol=1e-7,atol=1e-9)
  self.assertLess(np.linalg.norm(J-np.eye(3),axis=(1,2)).max(),field.lipschitz_bound)
  self.assertTrue(np.all(np.linalg.det(J)>0))
 def test_fixes_complete_triangle_and_retains_original_sides(self):
  t=np.array([[[-.02,-.02,0],[.02,-.02,0],[.02,.02,0]],[[-.02,-.02,0],[.02,.02,0],[-.02,.02,0]]])
  guard=TriangleSurfaceGuard.from_triangles(t,radius_m=.001)
  field=CompactTranslation([0,0,0],[.08,.08,.08],[0,0,-.012],'test')
  rng=np.random.default_rng(9);fixed=np.column_stack([rng.uniform(-.019,.019,(100,2)),np.zeros(100)])
  points=np.r_[fixed,[[.0011,.0022,.003],[.0011,.0022,-.003]]]
  q,J,p=guarded_translation_flow(points,[field],guard,fade_m=.003,steps=64)
  np.testing.assert_array_equal(q[:100],fixed)
  self.assertGreater(q[-2,2],0);self.assertLess(q[-2,2],points[-2,2]);self.assertLess(q[-1,2],0)
  self.assertTrue(np.all(np.linalg.det(J)>0));self.assertGreater(p['analytic_determinant_lower_bound'],0)
 def test_jacobian_matches_perturbation_away_from_cover_seams(self):
  t=np.array([[[0,0,0],[.01,0,0],[0,.01,0]]]);guard=TriangleSurfaceGuard.from_triangles(t,radius_m=.0005)
  field=CompactTranslation([0,0,.005],[.03,.03,.03],[.0002,0,-.0004],'test')
  x=np.array([[.00211,.00318,.00533]])
  q,J,_=guarded_translation_flow(x,[field],guard,fade_m=.002,steps=16)
  columns=[]
  for axis in range(3):
   e=np.eye(3)[axis]*1e-8
   a=guarded_translation_flow(x+e,[field],guard,fade_m=.002,steps=16)[0]
   b=guarded_translation_flow(x-e,[field],guard,fade_m=.002,steps=16)[0]
   columns.append(((a-b)/2e-8)[0])
  np.testing.assert_allclose(J[0],np.column_stack(columns),rtol=2e-6,atol=2e-7)
 def test_refuses_unbounded_cover_or_noninjective_step(self):
  t=np.array([[[0,0,0],[1,0,0],[0,1,0]]])
  with self.assertRaises(ValueError):TriangleSurfaceGuard.from_triangles(t,radius_m=.0001,maximum_samples=100)
  guard=TriangleSurfaceGuard(np.array([[0.,0,0]]),.001)
  f=CompactTranslation([0,0,0],[.1,.1,.1],[0,0,-.01],'test')
  with self.assertRaises(ValueError):guarded_translation_flow(np.zeros((1,3)),[f],guard,fade_m=.0001,steps=1)


class CorridorAndCropTests(unittest.TestCase):
 def test_corridor_permutation_and_gradient_bound(self):
  from geometry_tools.algorithms.surface_barrier_flow import SeparableCorridorField
  parameters=dict(centerline_stations=[dict(y_m=-.04,x_m=.02,z_m=.05),dict(y_m=.04,x_m=.02,z_m=.05)],center_y_m=0,core_axes_m=[.005,.02,.01],support_axes_m=[.02,.04,.03],translation_m=[-.0007,0,-.0005],y_core_interval_m=[-.02,.02],y_support_interval_m=[-.04,.04])
  field=SeparableCorridorField(parameters,[0,2,1]);q=np.array([[.02,.05,0],[.03,.06,.028]])
  mapped=field.map(q);J=field.jacobians(q)
  np.testing.assert_allclose(mapped[0]-q[0],[-.0007,-.0005,0],atol=1e-17)
  self.assertTrue(np.all(np.linalg.svd(J-np.eye(3),compute_uv=False)[:,0]<=field.lipschitz_bound))
  for axis in range(3):
   e=np.eye(3)[axis]*1e-8;actual=(field.map(q+e)-field.map(q-e))/2e-8
   np.testing.assert_allclose(actual,J[:,:,axis],atol=2e-8)
 def test_swept_crop_keeps_a_large_face_and_discards_only_unreachable_face(self):
  from geometry_tools.algorithms.surface_barrier_flow import swept_surface_guard_mask
  fixed=np.array([[[-1,-1,0],[1,-1,0],[0,1,0]],[[4,4,0],[5,4,0],[4,5,0]]],float)
  moving=np.array([[[-.01,-.01,.02],[.01,-.01,.02],[0,.01,.02]]])
  keep,proof=swept_surface_guard_mask(fixed,moving,[0,0,.03])
  np.testing.assert_array_equal(keep,[True,False]);self.assertTrue(proof['allPotentialMovingFaceContactsCovered'])


class AnisotropicCoverTests(unittest.TestCase):
 def test_long_thin_triangle_coverage_without_square_oversampling(self):
  t=np.array([[[0.,0,0],[.1,0,0],[.04,.00002,0]]])
  with self.assertRaises(ValueError):TriangleSurfaceGuard.from_triangles(t,radius_m=.001,maximum_samples=1000)
  guard=TriangleSurfaceGuard.from_triangles(t,radius_m=.001,maximum_samples=1000,anisotropic=True)
  rng=np.random.default_rng(19);weights=rng.dirichlet([1,1,1],size=2000);points=weights@t[0]
  distance,_=guard.lower_distance_and_gradient(points)
  np.testing.assert_array_equal(distance,np.zeros(len(points)));self.assertLess(len(guard.points),300);self.assertEqual(guard.evidence['anisotropicStripTriangles'],1)
 def test_oblique_and_degenerate_triangles_are_covered(self):
  triangles=np.array([[[.1,.2,.3],[.11,.21,.32],[.108,.207,.315]],[[0,0,0],[.1,0,0],[.06,0,0]],[[0,0,0],[0,0,0],[0,0,0]]])
  guard=TriangleSurfaceGuard.from_triangles(triangles,radius_m=.0005,maximum_samples=10000,anisotropic=True)
  rng=np.random.default_rng(3)
  for t in triangles:
   points=rng.dirichlet([1,1,1],size=1000)@t
   np.testing.assert_array_equal(guard.lower_distance_and_gradient(points)[0],np.zeros(1000))


class SkinBandTests(unittest.TestCase):
 def test_keeps_deep_controls_literal_and_transports_only_superficial_band(self):
  skin=np.array([[[-.03,-.03,0],[.03,-.03,0],[.03,.03,0]],[[-.03,-.03,0],[.03,.03,0],[-.03,.03,0]]])
  deep=skin.copy();deep[:,:,2]=-.02
  band=TriangleSurfaceGuard.from_triangles(skin,radius_m=.001,anisotropic=True)
  barrier=TriangleSurfaceGuard.from_triangles(deep,radius_m=.001,anisotropic=True)
  field=CompactTranslation([0,0,0],[.06,.06,.06],[0,0,-.001],'test')
  x=np.array([[.001,.002,-.012],[.001,.002,.0008]])
  q,J,p=guarded_translation_flow(x,[field],barrier,fade_m=.0015,steps=32,skin_band_guard=band,band_inner_m=.001,band_outer_m=.005)
  np.testing.assert_array_equal(q[0],x[0]);np.testing.assert_array_equal(J[0],np.eye(3));self.assertLess(q[1,2],x[1,2]);self.assertGreater(p['analytic_determinant_lower_bound'],0)
