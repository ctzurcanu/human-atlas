import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_collar_union import union_source_collars
from geometry_tools.algorithms.tube_continuity import continuity_audit
from geometry_tools.open3d_backend import diagnostics


class SourceCollarUnionTests(unittest.TestCase):
 def source(self):
  theta=np.arange(6)*np.pi/3;circle=np.c_[np.cos(theta),np.sin(theta),np.zeros(6)]*.00025;points=circle.tolist();faces=[]
  for direction in [np.array([0.,0.,-1.]),np.array([.5,0.,1.]),np.array([-.5,0.,1.])]:
   direction/=np.linalg.norm(direction);u=np.cross([0.,1.,0.],direction);u/=np.linalg.norm(u);v=np.cross(direction,u)
   if direction[2]<0:u=np.array([1.,0.,0.]);v=np.array([0.,1.,0.])
   previous=np.arange(6)
   for level in [1,2]:
    q=direction*.0015*level+np.outer(np.cos(theta),u)*.00025+np.outer(np.sin(theta),v)*.00025;ids=np.arange(len(points),len(points)+6);points.extend(q.tolist())
    for i,j in zip(range(6),np.roll(np.arange(6),-1)):
     pair=np.array([[previous[i],previous[j],ids[j]],[previous[i],ids[j],ids[i]]]);faces.extend((pair[:,[0,2,1]]if direction[2]<0 else pair).tolist())
    previous=ids
  return Mesh(np.array(points),np.array(faces),'analytical-metres')

 def test_three_shared_tubes_become_one_open_manifold_junction(self):
  m=self.source();q,p=union_source_collars(m);self.assertTrue(p['carrierAudit']['topological_continuity_pass']);self.assertEqual(p['computationRootCapsRetained'],0);self.assertEqual(p['emittedLumenDisks'],0);np.testing.assert_array_equal(q.vertices[:len(m.vertices)],m.vertices);self.assertTrue(continuity_audit(q.vertices,q.faces,expected_mouths=3,expected_components=1,expected_euler=-1)['topological_continuity_pass']);self.assertEqual(diagnostics(q)['self_intersection_pairs'],0)

 def test_unrelated_disconnected_source_tube_stays_literal(self):
  source=self.source();offset=len(source.vertices);vertices=np.r_[source.vertices,source.vertices[:18]+[.02,0,0]];extra=source.faces[:24]+offset;faces=np.r_[source.faces,extra];m=Mesh(vertices,faces,source.frame);q,p=union_source_collars(m);retained=np.array(p['retainedSourceFaceIds']);unrelated=np.flatnonzero(retained>=len(source.faces));np.testing.assert_array_equal(q.faces[unrelated],extra);np.testing.assert_array_equal(q.vertices[:len(vertices)],vertices);self.assertTrue(continuity_audit(q.vertices,q.faces,expected_mouths=5,expected_components=2,expected_euler=-1)['topological_continuity_pass']);self.assertEqual(diagnostics(q)['self_intersection_pairs'],0);self.assertFalse(set(p['removedSourceFaceIds'])&set(range(len(source.faces),len(faces))))

 def test_exact_position_seams_retain_distinct_uv_corners(self):
  source=self.source();vertices=source.triangles.reshape(-1,3);faces=np.arange(len(vertices)).reshape(-1,3);uvs=np.c_[np.arange(len(vertices)),np.arange(len(vertices))*.5].astype('<f4');m=Mesh(vertices,faces,source.frame,{'uvs':uvs});q,p=union_source_collars(m);np.testing.assert_array_equal(q.attributes['uvs'][:len(vertices)],uvs);retained=np.array(p['retainedSourceFaceIds']);np.testing.assert_array_equal(q.faces[:len(retained)],faces[retained]);self.assertTrue(continuity_audit(q.vertices,q.faces,expected_mouths=3,expected_components=1,expected_euler=-1)['topological_continuity_pass']);self.assertEqual(diagnostics(q)['self_intersection_pairs'],0)

if __name__=='__main__':unittest.main()
