import unittest
import numpy as np
from geometry_tools.algorithms.passage_cycles import handle_cycle


class PassageCycleTests(unittest.TestCase):
 def torus(self):
  n,m=40,12;u=np.arange(n)*2*np.pi/n;v=np.arange(m)*2*np.pi/m
  p=np.array([[(.015+.004*np.cos(b))*np.cos(a),.004*np.sin(b),(.015+.004*np.cos(b))*np.sin(a)]for a in u for b in v]);f=[]
  for i in range(n):
   for j in range(m):
    a=i*m+j;b=((i+1)%n)*m+j;c=((i+1)%n)*m+(j+1)%m;d=i*m+(j+1)%m;f.extend([[a,b,c],[a,c,d]])
  return p,np.asarray(f)
 def test_actual_handle_cycle_stays_on_existing_edges_and_is_nonseparating(self):
  p,f=self.torus();cycle,proof=handle_cycle(p,f,[0,0,0],radius=.025,half_length=.01)
  edges={tuple(sorted((a,b)))for face in f for a,b in zip(face,np.roll(face,-1))}
  self.assertTrue(all(tuple(sorted((a,b)))in edges for a,b in zip(cycle,np.roll(cycle,-1))))
  self.assertTrue(proof['nonseparating_source_cycle']);self.assertEqual(proof['face_dual_components_after_cut'],1)
  self.assertLess(proof['perimeter_m'],2*np.pi*.012)
 def test_wrong_conduit_axis_does_not_invent_a_surface_handle(self):
  p,f=self.torus()
  with self.assertRaisesRegex(ValueError,'No nonseparating'):handle_cycle(p,f,[.06,0,0],radius=.02,half_length=.01)
 def test_excluded_axis_cannot_be_hidden_inside_the_same_selected_cycle(self):
  p,f=self.torus()
  with self.assertRaisesRegex(ValueError,'No nonseparating'):handle_cycle(p,f,[0,0,0],radius=.025,half_length=.01,excluded_centers=[[0,0,0]])
 def test_winding_around_a_solid_sphere_is_not_mistaken_for_a_passage(self):
  p=.01*np.array([[1,0,0],[-1,0,0],[0,1,0],[0,-1,0],[0,0,1],[0,0,-1]],float)
  f=np.array([[2,0,4],[2,4,1],[2,1,5],[2,5,0],[3,4,0],[3,1,4],[3,5,1],[3,0,5]])
  with self.assertRaisesRegex(ValueError,'No nonseparating'):handle_cycle(p,f,[0,0,0],radius=.02,half_length=.02)

if __name__=='__main__':unittest.main()
