"""Close one declared source-surface defect without moving a leaflet control."""
from pathlib import Path
import argparse,json,sys,time,functools
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_npz,atomic_bytes,digest_file
from geometry_tools.audit import audit,topology
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from geometry_tools.algorithms.finite_cell_bank import pack_finite_cells
from geometry_tools.tissue import decode_packed
ap=argparse.ArgumentParser();ap.add_argument('--stage',required=True);args=ap.parse_args();start=time.monotonic();h=ROOT/'.local-models/geometry-workbench/heart';out=h/args.stage;out.mkdir(exist_ok=False);source=Mesh.load(h/'apparatus-source-v203/part-7.npz');reader=AtlasReader();current=reader.reload(source.source,source.frame);assert current.geometry_sha==source.geometry_sha
original=audit(source,require_closed=False,max_pairs=500000,max_tests=4000000);atomic_json(out/'source-audit.json',original);assert original['passed'];boundary=json.loads((h/'pulmonary-cusp-boundary-v222/boundary.json').read_text());p,inv=np.unique(source.vertices,axis=0,return_inverse=True);f=inv[source.faces];loop=tuple(boundary['boundary_vertex_ids'][::-1]);assert np.array_equal(p[np.array(loop[::-1])],np.array(boundary['boundary_m']));tree=BVH(source)
@functools.lru_cache(None)
def triangulations(ring):
 if len(ring)<3:return((),)
 if len(ring)==3:return((ring,),)
 result=[]
 for k in range(1,len(ring)-1):
  for left in triangulations(ring[:k+1]):
   for right in triangulations(ring[k:]):result.append((*left,*right,(ring[0],ring[k],ring[-1])))
 return tuple(result)
scores=[]
for tr in triangulations(loop):
 tri=np.array(tr,int);xyz=p[tri];cross=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);area=np.linalg.norm(cross,axis=1)*.5
 if np.any(area==0):continue
 scores.append((float(area.sum()),float(np.linalg.norm(np.diff(xyz[:,[0,1,2,0]],axis=1),axis=2).sum()),tri))
scores.sort(key=lambda row:row[:2]);attempts=[];chosen=None
for area,perimeter,patch in scores:
 if time.monotonic()-start>25:raise TimeoutError('Cusp patch search stage deadline')
 patchmesh=Mesh(p,patch,source.frame);_,contact=pair_contacts(patchmesh,source,0.,second_bvh=tree,max_pairs=50000,max_tests=500000);candidate=Mesh(p,np.concatenate([f,patch]),source.frame);t=topology(candidate)
 row=dict(area_m2=area,source_contact=contact,closed=t['closed_oriented_manifold'],genus=[g['closed_orientable_genus']for g in t['components']]);attempts.append(row)
 if not contact['complete_within_pair_budget']or contact['crossing_face_pairs']or contact['coplanar_overlap_face_pairs']or not t['closed_oriented_manifold']or row['genus']!=[0]:continue
 checked=audit(candidate,require_closed=True,expected_components=1,expected_genus=[0],max_pairs=500000,max_tests=4000000);row['full_candidate_audit']=checked
 if checked['passed']:chosen=patch;break
atomic_json(out/'search.json',dict(passed=chosen is not None,attempts=attempts,possible_triangulations=len(scores)));assert chosen is not None,'No finite patch passes the whole source checks'
record=source.source['record'];raw=reader._bank(Path(source.source['bank_path']));used=np.array(source.source['source_vertex_ids'],int);words=reader.view(raw,record['normals'],record['vertexCount'],'<i2',3)[used][source.faces].copy();xyz=p[f];cross=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0]);bad=np.einsum('tki,ti->tk',words.astype(float),cross)<=0;fixed=np.rint(cross/np.linalg.norm(cross,axis=1)[:,None]*32767).astype('<i2');words[bad]=np.broadcast_to(fixed[:,None],words.shape)[bad];patchtri=p[chosen];normal=np.cross(patchtri[:,1]-patchtri[:,0],patchtri[:,2]-patchtri[:,0]);normal/=np.linalg.norm(normal,axis=1)[:,None];patchwords=np.broadcast_to(np.rint(normal*32767).astype('<i2')[:,None],(len(chosen),3,3)).copy();uv=source.attributes['uvs'];uv_by_position=np.zeros((len(p),2),dtype='<f4');uv_by_position[inv]=uv;surfaces=dict(exterior=f,repair=chosen);fields=dict(exterior=words,repair=patchwords);uvs=dict(exterior=uv[source.faces],repair=uv_by_position[chosen]);pid=record['id'];packed=pack_finite_cells(p,surfaces,fields,{'cusp':[dict(surface='exterior',side=1),dict(surface='repair',side=1)]},{'cusp':pid},corner_uvs=uvs);atomic_bytes(out/'cells.bin',packed.data);decoded=decode_packed(packed,source.frame,np.zeros(3))['cusp'];checked=audit(decoded,require_closed=True,expected_components=1,expected_genus=[0],max_pairs=500000,max_tests=4000000);atomic_json(out/'checked-geometry.json',checked);assert checked['passed'];assert np.array_equal(decoded.triangles[:len(f)],source.triangles)
atomic_npz(out/'source-arrangement.npz',positions=p,source_faces=f,patch=chosen,source_normal_words=words,patch_normal_words=patchwords,source_corner_uvs=uv[source.faces],patch_corner_uvs=uv_by_position[chosen]);emission=dict(passed=True,source_part_id=pid,source_geometry_sha256=source.geometry_sha,source_attributes_sha256=source.attributes_sha,source_record=record,parts=packed.cells,bank_sha256=digest_file(out/'cells.bin'),bytes=len(packed.data),source_face_count=len(f),patch_triangles=len(chosen),source_normal_corners_repaired=int(bad.sum()),original_source_triangles_literal=True,original_uvs_literal=True,new_points=0,patch_area_m2=attempts[-1]['area_m2'],elapsed_seconds=time.monotonic()-start,installed=False,anatomical_acceptance=False);atomic_json(out/'emission.json',emission);print(json.dumps({k:v for k,v in emission.items()if k not in ['source_record','parts']}),flush=True)
