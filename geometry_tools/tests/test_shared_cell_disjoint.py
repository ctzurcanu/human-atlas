import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.shared_cell_disjoint import shared_cell_disjoint

def tetra(points):
 p=np.asarray(points,float);f=np.array([[0,2,1],[0,1,3],[1,2,3],[2,0,3]])
 if np.linalg.det((p[1:]-p[0]).T)<0:f=f[:,[0,2,1]]
 return Mesh(p,f,'test')

class SharedCellDisjointTests(unittest.TestCase):
 def setUp(self):self.a=tetra([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
 def test_opposed_common_facet(self):
  b=tetra([[0,0,0],[1,0,0],[0,1,0],[0,0,-1]])
  self.assertTrue(shared_cell_disjoint(self.a,b)['passed'])
 def test_identical_occupied_cells_rejected(self):
  self.assertFalse(shared_cell_disjoint(self.a,self.a)['passed'])
 def test_contained_cell_rejected(self):
  b=tetra([[.1,.1,.1],[.2,.1,.1],[.1,.2,.1],[.1,.1,.2]])
  self.assertFalse(shared_cell_disjoint(self.a,b)['passed'])
 def test_crossing_cells_rejected(self):
  b=tetra([[.2,.2,-.2],[1.2,.2,-.2],[.2,1.2,-.2],[.2,.2,.8]])
  self.assertFalse(shared_cell_disjoint(self.a,b)['passed'])
 def test_adjacent_component_does_not_hide_an_embedded_component(self):
  outside=tetra([[0,0,0],[1,0,0],[0,1,0],[0,0,-1]])
  inside=tetra([[.1,.1,.1],[.2,.1,.1],[.1,.2,.1],[.1,.1,.2]])
  b=Mesh(np.r_[outside.vertices,inside.vertices],np.r_[outside.faces,inside.faces+4],'test')
  self.assertFalse(shared_cell_disjoint(self.a,b)['passed'])

if __name__=='__main__':unittest.main()
