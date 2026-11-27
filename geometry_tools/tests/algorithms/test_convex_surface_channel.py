"""Analytical evidence for a shared conforming polygonal opening.

Run with the task Python runtime. These fixtures test hole topology, projected
area and source-control preservation, independently of thyroid shape choices.
"""
import unittest
import numpy as np
from geometry_tools.algorithms.convex_surface_channel import PointBank,subtract_surface,area
from scripts.ta98_geometry_checks import topology

class ConvexSurfaceChannelTests(unittest.TestCase):
 def setUp(self):
  self.points=np.array([[x,y,0] for y in np.linspace(0,1,7) for x in np.linspace(0,1,7)],np.float32)
  faces=[]
  for y in range(6):
   for x in range(6):
    a=y*7+x;faces.extend([[a,a+1,a+8],[a,a+8,a+7]])
  self.faces=np.array(faces,np.uint32);self.normals=np.tile(np.array([0,0,1],np.float32),(len(self.points),1))

 def test_rotated_hole_preserves_shared_source_plane_and_area(self):
  angle=.57;rotation=np.array([[np.cos(angle),-np.sin(angle)],[np.sin(angle),np.cos(angle)]])
  polygon=np.array([[-.17,-.11],[.17,-.11],[.17,.11],[-.17,.11]])@rotation.T
  bank=PointBank(self.points,self.normals);faces,parents,detail=subtract_surface(bank,self.faces,np.array([.5,.5,0]),np.eye(3),polygon,(-.01,.01),'fixture')
  positions=np.array(bank.positions,np.float32);t=topology(positions,faces);tri=positions[faces].astype(float);actual=np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1).sum()/2
  self.assertEqual(len(t['components']),1);self.assertEqual(len(t['boundaryLoops']),2)
  for key in ['nonmanifoldEdges','inconsistentWindingEdges','duplicateFaces','degenerateFaces']:self.assertEqual(t[key],0)
  self.assertTrue(all(q['cycle'] for q in t['boundaryLoops']));self.assertEqual(t['invalidCycleOrPathVertexLinks'],[])
  self.assertAlmostEqual(actual,1-abs(area(polygon)),places=7)
  self.assertTrue(np.array_equal(positions[:len(self.points)],self.points));self.assertTrue(np.array_equal(np.array(bank.normals[:len(self.points)],np.float32),self.normals))
  self.assertTrue(np.all(positions[:,2]==0));self.assertTrue(len(detail['sourceFacesCut'])>0)
  for parent in set(range(len(self.faces)))-set(detail['sourceFacesCut'])-set(detail['adjacentSourceFacesRetessellated']):
   self.assertTrue(any(np.array_equal(tri,self.faces[parent]) for tri in faces[parents==parent]))

 def test_missing_cutter_keeps_every_face_and_control_literal(self):
  bank=PointBank(self.points,self.normals);polygon=np.array([[-.1,-.1],[.1,-.1],[.1,.1],[-.1,.1]])
  faces,parents,detail=subtract_surface(bank,self.faces,np.array([2,2,0]),np.eye(3),polygon,(-.01,.01),'miss')
  self.assertTrue(np.array_equal(faces,self.faces));self.assertTrue(np.array_equal(parents,np.arange(len(self.faces))))
  self.assertEqual(len(bank.positions),len(self.points));self.assertEqual(detail['sourceFacesCut'],[])

if __name__=='__main__':unittest.main()
