"""Actual section boundaries must not turn into bounding boxes or fake caps."""
import importlib.util
from pathlib import Path
import unittest
import numpy as np

spec=importlib.util.spec_from_file_location('tibial_corridor_sections',Path(__file__).resolve().parents[2]/'scripts/derive-tibial-compartment-course-v39.py')
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)

def cut(points,opened=False):
 p=np.asarray(points,float);v=np.c_[p[:,0],np.zeros(len(p)),p[:,1]];segments=np.stack([v,np.roll(v,-1,axis=0)],axis=1)
 return dict(segments=segments[:-1]if opened else segments,loops=[]if opened else[v],open_components=[{}]if opened else[])

class ActualCorridorSections(unittest.TestCase):
 def test_concavity_remains_available(self):
  c=cut([[0,0],[2,0],[2,1],[1,1],[1,2],[0,2]])
  self.assertEqual(module.outside_bone(c,[[1.5,1.5],[.5,.5]],.001).tolist(),[True,False])

 def test_real_boundary_margin(self):
  c=cut([[0,0],[1,0],[1,1],[0,1]])
  self.assertEqual(module.outside_bone(c,[[1.0002,.5],[1.002,.5],[1,.5]],.001).tolist(),[False,True,False])

 def test_open_contour_is_not_silently_closed(self):
  c=cut([[0,0],[1,0],[1,1],[0,1]],opened=True)
  self.assertEqual(module.outside_bone(c,[[-1,.5],[.5,-1],[1.1,.5],[.5,1.1]],.001).tolist(),[False,False,True,True])

 def test_missing_fascial_chart_stays_missing(self):
  segments=np.array([[[0,0,1],[1,0,1]],[[2,0,1],[3,0,1]]],float)
  result=module.front_cover(segments,np.array([.5,1.5,2.5]))
  self.assertEqual(result[[0,2]].tolist(),[1,1]);self.assertTrue(np.isneginf(result[1]))

if __name__=='__main__':unittest.main()
