import unittest
import numpy as np
from geometry_tools.algorithms.source_root_sections import recover_root_sections


class SourceRootSectionsTest(unittest.TestCase):
 def test_twelve_to_six_neck_retains_literal_faces_without_exporter_quads(self):
  p=[]
  for y,n in [(0,12),(1,12),(2,6),(3,6)]:
   p.extend([[np.cos(t),y,np.sin(t)]for t in np.arange(n)*2*np.pi/n])
  f=[]
  for i in range(12):
   j=(i+1)%12;f.extend([[i,12+i,12+j],[i,12+j,j]])
  for i in range(6):
   a=12+2*i;b=12+(2*i+1)%12;c=12+(2*i+2)%12;d=24+i;e=24+(i+1)%6;f.extend([[a,d,b],[b,d,e],[b,e,c]])
  for i in range(6):
   j=(i+1)%6;f.extend([[24+i,30+i,30+j],[24+i,30+j,24+j]])
  p=np.array(p,float);f=np.array(f,int);oldp=p.copy();oldf=f.copy();r=recover_root_sections(p,f,np.arange(12))
  self.assertEqual([len(x['rawIds'])for x in r],[12,6,6]);self.assertTrue(all(x['certifiedManifoldAnnulus']for x in r))
  self.assertEqual(sum(len(x['sourceFaceIds'])for x in r),len(f));np.testing.assert_array_equal(p,oldp);np.testing.assert_array_equal(f,oldf)

if __name__=='__main__':unittest.main()
