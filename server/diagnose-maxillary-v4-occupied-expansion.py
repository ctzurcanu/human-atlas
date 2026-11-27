"""Bounded source-classifier occupancy evidence; not an anatomical acceptance."""
from pathlib import Path
import json,time,sys,hashlib
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology,winding_contains
from geometry_tools.algorithms.annular_pipe_envelopes import boundary_loops
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/maxillary-native-v1/junction-repair-v1';C=B/'candidate-v4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def closed_course(sm,rr):
 ids=np.flatnonzero(np.all(np.isin(sm.faces,rr.ravel()),axis=1));faces=sm.faces[ids];verts=sm.vertices.tolist();out=faces.tolist();caps=[]
 for loop in boundary_loops(faces):
  centre=len(verts);verts.append(sm.vertices[loop].mean(0).tolist())
  for a,b in zip(loop,np.roll(loop,-1)):caps.append([int(b),int(a),centre]);out.append(caps[-1])
 m=Mesh(np.array(verts),np.array(out),sm.frame,{},dict(diagnostic_only=True));assert topology(m)['closed_oriented_manifold'];return m,Mesh(m.vertices,np.array(caps),m.frame,{},m.source)
def main():
 start=time.monotonic();p=json.loads((B/'plan.json').read_text());contract=json.loads((B/'local-completion-contract-v4.json').read_text());cp=ROOT/p['source_catalogue'];assert sha(cp)==p['source_master_sha256'];cat=json.loads(cp.read_text());src={r['source']:r for r in p['sources']};a={k:v for k,v in np.load(B/'source-course-arrays.npz').items()};an={k:v for k,v in np.load(C/'ancestry.npz').items()};m=Mesh.load(C/'cohort.npz');reader=AtlasReader();rows=[]
 for sid,side in[(0,'l'),(1,'r')]:
  sm=Mesh.load(src[sid]['path']);c6,cap6=closed_course(sm,a[f's{sid}_c6_rings']);c4,cap4=closed_course(sm,a[f's{sid}_c4_rings']);anchor=np.array(next(x['anchor_m']for x in contract['anchors']if x['source']==sid and x['component']==6));patchids=[i for i in np.flatnonzero((an['source_indices']==sid)&(an['source_face_ids']==-2))if any(np.array_equal(v,anchor)for v in an['parent_triangles'][i])];tri=m.triangles[patchids]
  # Uniform subdivision gives a reproducible Lipschitz upper bound over whole triangles.
  tiny=tri.copy()
  for level in range(3):
   u,v,w=tiny[:,0],tiny[:,1],tiny[:,2];uv=(u+v)/2;vw=(v+w)/2;wu=(w+u)/2;tiny=np.concatenate([np.stack([u,uv,wu],1),np.stack([uv,v,vw],1),np.stack([wu,vw,w],1),np.stack([uv,vw,wu],1)])
  pts=tiny.mean(1);cover=np.linalg.norm(tiny-pts[:,None,:],axis=2).max(1);inside=[];dist=[]
  for old in[c6,c4]:
   evidence=winding_contains(old,pts,boundary_tolerance_m=1e-12);inside.append(np.array(evidence['inside'],bool)|np.array(evidence['boundary_ambiguous'],bool));nearest,_=BVH(old).closest_points(pts);dist.append(nearest['distance_m'])
  old_union_inside=inside[0]|inside[1];outside=np.where(old_union_inside,0,np.minimum(dist[0],dist[1]));q=int(np.argmax(outside));tooth_id=f'ZA:Upper second molar tooth.{side}:Teeth-roots';tooth=reader.read_part(cat,next(x for x in cat['parts']if x['id']==tooth_id),frame=m.frame,catalogue_path=cp,catalogue_sha=p['source_master_sha256']);caprows=[]
  for component,caps in[(6,cap6),(4,cap4)]:
   pairs,summary=pair_contacts(caps,tooth,0.,max_pairs=20000,max_tests=100000);assert summary['complete_within_pair_budget'];bad=pairs['pairs'][np.isin(pairs['relation'],['crossing','coplanar-overlap'])];caprows.append(dict(component=component,strict_pairs=len(bad),tooth_faces=sorted(set(map(int,bad[:,1])))))
  rows.append(dict(source=sid,candidate_new_c6_sleeve_faces=patchids,triangle_count=len(tri),subdivision_levels=3,whole_triangle_subcells=len(tiny),outside_probe_count=int((outside>0).sum()),maximum_sampled_outside_m=float(outside.max()),conservative_whole_triangle_outside_upper_bound_m=float(np.max(outside+cover)),maximum_cover_radius_m=float(cover.max()),maximum_outside_witness_m=pts[q].tolist(),old_capped_classifiers=[dict(component=6,geometry_sha256=c6.geometry_sha),dict(component=4,geometry_sha256=c4.geometry_sha)],old_cap_tooth_contacts=caprows,interpretation='Outside distance is to the union of the two original closed source classifiers. Surface samples prove occupied expansion when positive; cover radius gives a conservative whole-triangle upper bound. Neither inherited nor new tooth contact is accepted.'))
 report=dict(schema='maxillary-v4-source-classifier-occupied-expansion/v1',complete=True,source_master_sha256=p['source_master_sha256'],candidate_sha256=sha(C/'cohort.npz'),rows=rows,elapsed_seconds=time.monotonic()-start,anatomical_acceptance=False,geometry_modified=False);serial=lambda x:x.tolist()if hasattr(x,'tolist')else x.item();(C/'source-classifier-occupied-expansion.json').write_text(json.dumps(report,indent=2,default=serial));print(json.dumps(report,default=serial))
if __name__=='__main__':main()
