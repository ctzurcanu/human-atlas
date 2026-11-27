import json,tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.native_labelled_tetrahedra import partition_constraints_native
from geometry_tools.algorithms.shared_float32_boundary import normalize_shared_boundary
from geometry_tools.model import Mesh
from geometry_tools.audit import audit

class SharedFloat32BoundaryTests(unittest.TestCase):
 def test_coincident_owned_plane_fixture_is_closed_after_actual_encoding(self):
  fixture=json.loads((Path(__file__).with_name('fixtures')/'coincident_owned_planes.json').read_text());p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]])
  with tempfile.TemporaryDirectory()as folder:r=partition_constraints_native(p,np.array([[0,1,2,3]]),fixture['fields'],fixture['cells'],fixture['zero_side'],Path(folder)/'native',seconds=5)
  encoded=normalize_shared_boundary(r['positions'],r['cells']);total=0.;owners={}
  for name,faces in encoded['cells'].items():
   if not len(faces):continue
   result=audit(Mesh(encoded['positions'],faces,'test'),require_closed=True);self.assertTrue(result['passed'],result['failures']);total+=result['topology']['oriented_triangle_integral_m3']
   for f in faces:owners.setdefault(tuple(sorted(map(int,f))),[]).append(name)
  self.assertAlmostEqual(total,1/6,places=7);self.assertTrue(all(len(x)<=2 for x in owners.values()))
 def test_duplicate_outward_tissue_faces_are_rejected(self):
  p=np.eye(3)
  with self.assertRaisesRegex(ValueError,'Identically oriented'):normalize_shared_boundary(p,{'bad':np.array([[0,1,2],[0,1,2]])})
 def test_reciprocal_cancellation_is_independent_of_facet_order(self):
  p=np.eye(3);positive=[0,1,2];negative=[0,2,1]
  for faces in [[positive,positive,negative],[negative,positive,positive],[positive,negative,positive]]:
   result=normalize_shared_boundary(p,{'cell':np.array(faces)})
   self.assertEqual(len(result['cells']['cell']),1)
   self.assertEqual(result['evidence']['cells']['cell']['reciprocal_calculation_triangles_cancelled'],2)
   q=result['positions'][result['cells']['cell'][0]]
   self.assertGreater(np.dot(np.cross(q[1]-q[0],q[2]-q[0]),np.cross(p[1]-p[0],p[2]-p[0])),0)
