import unittest
import numpy as np
from geometry_tools.algorithms.source_annular_refinement import refine_annuli
from geometry_tools.algorithms.source_tube_topology import recover_source_annuli
from geometry_tools.model import Mesh
from geometry_tools.audit import topology
from geometry_tools.spatial import BVH
from geometry_tools.distance import triangle_pairs

def tube(stations,radius=.045):
 angle=np.arange(6)*2*np.pi/6;ring=np.c_[radius*np.cos(angle),np.zeros(6),radius*np.sin(angle)];v=np.concatenate([ring+[0,t,0]for t in stations]);faces=[]
 for k in range(len(stations)-1):
  for i in range(6):a=k*6+i;b=k*6+(i+1)%6;c=b+6;d=a+6;faces.extend([[a,b,c],[a,c,d]])
 return v,np.asarray(faces)

def curve(points,radius=1):
 t=points[:,1];return np.c_[radius*np.cos(t*np.pi/2),np.zeros(len(t)),radius*np.sin(t*np.pi/2)]

def crosses(v,f,w,g):
 a,b=Mesh(v,f,'fixture-metres'),Mesh(w,g,'fixture-metres');at,bt=a.triangles,b.triangles
 for pairs in BVH(a).candidate_batches(BVH(b),0,batch_size=1024):
  if np.any(np.isin(triangle_pairs(at[pairs[:,0]],bt[pairs[:,1]])['relation'],['crossing','coplanar-overlap'])):return True
 return False

class SourceAnnularRefinementTest(unittest.TestCase):
 def candidate(self,extra=False):
  v,f=tube([0,1]);pair=np.arange(len(f)).reshape(-1,2)
  if extra:v=np.r_[v,[[4,0,0],[4,1,0],[4,0,1]]];f=np.r_[f,[[12,13,14]]]
  ann=recover_source_annuli(v,f,pair);q=v+curve(v);result=refine_annuli(v,q,f,ann,pair,[0],curve,maximum_station_m=.01);return v,q,f,result
 def test_original_profiles_remain_and_source_annulus_stays_open(self):
  v,q,f,r=self.candidate();np.testing.assert_array_equal(r['vertices'][:len(q)],q);np.testing.assert_array_equal(r['source_vertices'][:len(v)],v)
  before=topology(Mesh(v,f,'fixture-metres'));after=topology(Mesh(r['vertices'],r['faces'],'fixture-metres'));self.assertEqual(after['open_edges'],before['open_edges']);self.assertEqual(after['nonmanifold_edges'],0);self.assertEqual(after['inconsistent_winding_edges'],0);self.assertEqual(sum(c['euler']for c in after['components']),sum(c['euler']for c in before['components']))
  for profile in r['new_profiles']:
   ids=profile['raw_ids'];a,b=r['source_vertices'][ids],r['vertices'][ids];np.testing.assert_allclose(np.linalg.norm(a[:,None]-a[None,:],axis=2),np.linalg.norm(b[:,None]-b[None,:],axis=2),atol=3e-16)
 def test_curved_station_surface_removes_a_coarse_chord_crossing(self):
  v,q,f,r=self.candidate();w,g=tube([0,.5,1]);w+=curve(w,.8);self.assertTrue(crosses(q,f,w,g));self.assertFalse(crosses(r['vertices'],r['faces'],w,g))
 def test_unassigned_source_faces_are_retained_with_one_owner(self):
  v,q,f,r=self.candidate(True);self.assertIn(12,r['kept_source_face_ids']);self.assertTrue(np.any(np.all(r['faces']==f[-1],axis=1)));self.assertEqual(len(r['replaced_source_face_ids']),12);self.assertEqual(len(set(r['replaced_source_face_ids'])),12)
 def test_duplicate_strip_and_vertex_budget_are_rejected(self):
  v,f=tube([0,1]);ann=recover_source_annuli(v,f,np.arange(len(f)).reshape(-1,2));q=v+curve(v)
  with self.assertRaises(ValueError):refine_annuli(v,q,f,ann,np.arange(len(f)).reshape(-1,2),[0,0],curve)
  with self.assertRaisesRegex(ValueError,'budget'):refine_annuli(v,q,f,ann,np.arange(len(f)).reshape(-1,2),[0],curve,maximum_new_vertices=5)
 def test_rotated_end_profiles_use_proper_intermediate_poses(self):
  from scipy.spatial.transform import Rotation
  v,f=tube([0,.02],radius=.001);ann=recover_source_annuli(v,f,np.arange(len(f)).reshape(-1,2));q=v.copy();center=v[6:].mean(0);q[6:]=(v[6:]-center)@Rotation.from_euler('x',45,degrees=True).as_matrix().T+center
  r=refine_annuli(v,q,f,ann,np.arange(len(f)).reshape(-1,2),[0],lambda p:np.zeros_like(p),maximum_station_m=.002,proper_profile_poses=True)
  np.testing.assert_array_equal(r['vertices'][:len(q)],q)
  for profile in r['new_profiles']:
   ids=profile['raw_ids'];a,b=r['source_vertices'][ids],r['vertices'][ids];np.testing.assert_allclose(np.linalg.norm(a[:,None]-a[None,:],axis=2),np.linalg.norm(b[:,None]-b[None,:],axis=2),atol=1e-16)
 def test_source_frame_blending_cannot_silently_shrink_calibre(self):
  from scipy.spatial.transform import Rotation
  v,f=tube([0,.02],radius=.001);center=v[6:].mean(0);v[6:]=(v[6:]-center)@Rotation.from_euler('x',90,degrees=True).as_matrix().T+center;ann=recover_source_annuli(v,f,np.arange(len(f)).reshape(-1,2))
  with self.assertRaisesRegex(ValueError,'area'):
   refine_annuli(v,v.copy(),f,ann,np.arange(len(f)).reshape(-1,2),[0],lambda p:np.zeros_like(p),maximum_station_m=.002,proper_profile_poses=True)
 def test_station_roll_is_retained_even_when_both_end_poses_are_identity(self):
  from scipy.spatial.transform import Rotation
  v,f=tube([0,.02],radius=.001);pairs=np.arange(len(f)).reshape(-1,2);ann=recover_source_annuli(v,f,pairs)
  def pose(points):return Rotation.from_euler('y',(45*np.sin(np.pi*points[:,1]/.02)**2)[:,None],degrees=True).as_matrix(),np.zeros_like(points)
  r=refine_annuli(v,v.copy(),f,ann,pairs,[0],lambda p:np.zeros_like(p),maximum_station_m=.002,proper_profile_poses=True,profile_pose=pose);np.testing.assert_array_equal(r['vertices'][:len(v)],v);middle=next(p for p in r['new_profiles']if abs(p['t']-.5)<1e-12);ids=middle['raw_ids'];expected=r['source_vertices'][ids]@Rotation.from_euler('y',45,degrees=True).as_matrix().T;np.testing.assert_allclose(r['vertices'][ids],expected,atol=1e-16)
 def test_already_short_explicit_branch_strips_retain_all_source_faces(self):
  v,f=tube([0,.001],radius=.0001);pairs=np.arange(len(f)).reshape(-1,2);ann=recover_source_annuli(v,f,pairs);r=refine_annuli(v,v.copy(),f,ann,pairs,[0],lambda p:np.zeros_like(p),maximum_station_m=.002,skip_short_strips=True);np.testing.assert_array_equal(r['faces'],f);np.testing.assert_array_equal(r['vertices'],v);self.assertEqual(r['new_triangles'],0);self.assertEqual(len(r['replaced_source_face_ids']),0)
 def test_proper_source_frames_retain_calibre_across_a_sharp_turn(self):
  from scipy.spatial.transform import Rotation
  v,f=tube([0,.02],radius=.001);centre=v[6:].mean(0);v[6:]=(v[6:]-centre)@Rotation.from_euler('x',90,degrees=True).as_matrix().T+centre;pairs=np.arange(len(f)).reshape(-1,2);ann=recover_source_annuli(v,f,pairs);r=refine_annuli(v,v.copy(),f,ann,pairs,[0],lambda p:np.zeros_like(p),maximum_station_m=.002,proper_profile_poses=True,proper_source_profile_strips=[0]);np.testing.assert_array_equal(r['vertices'][:len(v)],v)
  for profile in r['new_profiles']:
   p=r['source_vertices'][profile['raw_ids']];p=p-p.mean(0);area=.5*np.linalg.norm(np.cross(p,np.roll(p,-1,axis=0)).sum(0));self.assertAlmostEqual(area,3*np.sqrt(3)/2*.001**2,places=15);self.assertTrue(profile['properSourceProfilePose'])

if __name__=='__main__':unittest.main()
