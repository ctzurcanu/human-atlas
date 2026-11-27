import tempfile,unittest
from pathlib import Path
import numpy as np
from geometry_tools.native_labelled_tetrahedra import partition_native,partition_constraints_native
from geometry_tools.algorithms.labelled_tetrahedra import partition_labelled_tetrahedra
from geometry_tools.model import Mesh
from geometry_tools.audit import topology

class NativePartitionTests(unittest.TestCase):
 def fixture(self):
  p=np.array([[0,0,0],[1,0,0],[1,1,0],[0,1,0],[0,0,1],[1,0,1],[1,1,1],[0,1,1]],float)
  t=np.array([[0,1,2,6],[0,2,3,6],[0,3,7,6],[0,7,4,6],[0,4,5,6],[0,5,1,6]])
  return p,t
 def compare(self,scores):
  p,t=self.fixture();reference=partition_labelled_tetrahedra(p,t,scores,identical_score_priority=list(range(scores.shape[1])))
  with tempfile.TemporaryDirectory()as folder:result=partition_native(p,t,scores,Path(folder)/'native',seconds=5)
  volume=0
  for label,faces in result['cells'].items():
   rf=reference['cells'][label]
   if not len(faces):self.assertEqual(len(rf),0);continue
   m=Mesh(result['positions'],faces,'test');r=Mesh(reference['positions'],rf,'test');a=topology(m);b=topology(r)
   self.assertTrue(a['closed_oriented_manifold']);self.assertAlmostEqual(a['oriented_triangle_integral_m3'],b['oriented_triangle_integral_m3'],places=12);volume+=a['oriented_triangle_integral_m3']
  self.assertAlmostEqual(volume,1.,places=12)
 def test_two_regions(self):
  p,_=self.fixture();self.compare(np.c_[p[:,0]-.4,.4-p[:,0]])
 def test_three_regions_and_shared_cuts(self):
  p,_=self.fixture();self.compare(np.c_[p[:,0]-.3,p[:,1]-.4,p[:,2]-.2])
 def test_identical_affine_priority(self):
  p,_=self.fixture();self.compare(np.ones((len(p),2)))
 def test_layers_and_regional_core_share_interfaces(self):
  p,t=self.fixture();fields=np.c_[.1-p[:,0],p[:,0]-.9,p[:,1]-.5]
  cells=[[(0,1)],[(1,1)],[(0,-1),(1,-1),(2,-1)],[(0,-1),(1,-1),(2,1)]]
  with tempfile.TemporaryDirectory()as folder:r=partition_constraints_native(p,t,fields,cells,[-1,-1,1],Path(folder)/'layers',seconds=5)
  volumes=[];uses={}
  for owner,faces in r['cells'].items():
   proof=topology(Mesh(r['positions'],faces,'test'));self.assertTrue(proof['closed_oriented_manifold']);volumes.append(proof['oriented_triangle_integral_m3'])
   for face in faces:uses.setdefault(tuple(sorted(map(int,face))),[]).append(owner)
  np.testing.assert_allclose(volumes,[.1,.1,.4,.4],atol=1e-12);self.assertTrue(any(len(x)==2 for x in uses.values()));self.assertTrue(all(len(x)<=2 for x in uses.values()))
 def test_redundant_coincident_planes(self):
  p,t=self.fixture();a=.37-p[:,0]-2*p[:,1];fields=np.c_[a,2*a,p[:,2]-.41];cells=[[(0,1),(1,1),(2,1)],[(0,1),(1,1),(2,-1)],[(0,-1),(1,-1)]]
  with tempfile.TemporaryDirectory()as folder:r=partition_constraints_native(p,t,fields,cells,[1,1,1],Path(folder)/'redundant',seconds=5)
  total=0
  for faces in r['cells'].values():
   a=topology(Mesh(r['positions'],faces,'test'));self.assertTrue(a['closed_oriented_manifold']);total+=a['oriented_triangle_integral_m3']
  self.assertAlmostEqual(total,1.,places=12)
 def test_material_domains_keep_one_common_cut_at_their_interface(self):
  p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1],[0,0,-1]])
  t=np.array([[0,1,2,3],[0,2,1,4]]);fields=(p[:,0]-.25)[:,None];cells=[[(0,-1)],[(0,1)],[(0,-1)],[(0,1)]]
  with tempfile.TemporaryDirectory()as folder:r=partition_constraints_native(p,t,fields,cells,[1],Path(folder)/'domains',seconds=5,tetrahedron_domains=np.array([0,1]),cell_domains=np.array([0,0,1,1]))
  total=0.;uses={}
  for owner,faces in r['cells'].items():
   proof=topology(Mesh(r['positions'],faces,'test'));self.assertTrue(proof['closed_oriented_manifold']);total+=proof['oriented_triangle_integral_m3']
   for face in faces:uses.setdefault(tuple(sorted(map(int,face))),[]).append(owner)
  self.assertAlmostEqual(total,1/3,places=12);self.assertTrue(any(set(owners)=={0,2}for owners in uses.values()));self.assertTrue(any(set(owners)=={1,3}for owners in uses.values()));self.assertTrue(all(len(owners)<=2 for owners in uses.values()))
 def test_fields_coincident_only_on_a_source_face_keep_one_cut(self):
  p=np.array([[0.,0,0],[1,0,0],[0,1,0],[0,0,1]]);t=np.array([[0,1,2,3]]);fields=np.array([[0,.001291,.001291],[0,.001290,.001290],[0,-.000008,-.000008],[1,.00002,-.00003]]);cells=[[(1,a),(2,b)]for a,b in [(1,1),(1,-1),(-1,1),(-1,-1)]]
  with tempfile.TemporaryDirectory()as folder:r=partition_constraints_native(p,t,fields,cells,[1,1,1],Path(folder)/'face-coincidence',seconds=5)
  total=0.
  for faces in r['cells'].values():
   if not len(faces):continue
   proof=topology(Mesh(r['positions'],faces,'test'));self.assertTrue(proof['closed_oriented_manifold']);total+=proof['oriented_triangle_integral_m3']
  self.assertAlmostEqual(total,1/6,places=12)
 def test_nested_subdivision_refines_both_uses_of_shared_surface(self):
  p,t=self.fixture();fields=np.c_[.5-p[:,0],p[:,2]-.5];cells=[[(0,1),(1,1)],[(0,1),(1,-1)],[(0,-1)]]
  with tempfile.TemporaryDirectory()as folder:r=partition_constraints_native(p,t,fields,cells,[1,1],Path(folder)/'nested',seconds=5,surface_refinement_planes=[1])
  owners={}
  for owner,faces in r['cells'].items():
   self.assertTrue(topology(Mesh(r['positions'],faces,'test'))['closed_oriented_manifold'])
   for face in faces:owners.setdefault(tuple(sorted(map(int,face))),[]).append(owner)
  for key,uses in owners.items():
   q=r['positions'][list(key)];exterior=any(np.all(q[:,axis]==limit)for axis in range(3)for limit in [0,1]);self.assertEqual(len(uses),1 if exterior else 2)

if __name__=='__main__':unittest.main()
