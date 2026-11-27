import unittest
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.algorithms.source_stem_loft import replace_source_stem
from geometry_tools.algorithms.source_profile_axis import source_profile_graph
from geometry_tools.algorithms.port_hermite_course import cubic_metrics
from geometry_tools.algorithms.tube_continuity import continuity_audit


class SourceStemLoftTest(unittest.TestCase):
 def test_literal_mouths_retained_distal_faces_and_open_channel(self):
  theta=np.arange(6)*np.pi/3;p=np.array([[.002*np.cos(t),y,.002*np.sin(t)]for y in np.linspace(0,.05,7)for t in theta]);f=[]
  for k in range(6):
   for i in range(6):
    j=(i+1)%6;a,b,c,d=k*6+i,k*6+j,(k+1)*6+j,(k+1)*6+i;f.extend([[a,b,c],[a,c,d]])
  f=np.array(f,int);m=Mesh(p,f,'test:metres');g=source_profile_graph(p,f,triangle_pairing=np.arange(len(f)).reshape(-1,2));root=int(np.argmin([p[r,1].mean()for r in g['boundary_raw']]));T=np.diag([1.1,1.1,1.1,1.]);T[:3,3]=[.004,.01,.003]
  q,proof=replace_source_stem(m,root,np.eye(4),T)
  ids=proof['sourceRootRawIds'];np.testing.assert_array_equal(q.vertices[ids],p[ids]);retained=proof['retainedSourceFaceIds'];np.testing.assert_array_equal(q.triangles[:len(retained)],(p@T[:3,:3].T+T[:3,3])[f[retained]])
  audit=continuity_audit(q.vertices,q.faces,expected_mouths=2,expected_components=1,expected_euler=0);self.assertTrue(audit['topological_continuity_pass']);self.assertEqual(proof['lumenDisksCreated'],0);np.testing.assert_array_equal(m.vertices,p)
  areas=proof['sectionAreas_m2'];self.assertGreater(min(areas),0);self.assertTrue(all(min(areas[0],areas[-1])*.999<=x<=max(areas[0],areas[-1])*1.001 for x in areas))

 def test_straight_course_polynomial_handles_zero_components(self):
  metrics=cubic_metrics(np.array([[0.,0,0],[0,.04,0],[0,0,0],[0,0,0]]));self.assertGreater(metrics['minimumRadius_m'],1e20);self.assertAlmostEqual(metrics['length_m'],.04)

if __name__=='__main__':unittest.main()
