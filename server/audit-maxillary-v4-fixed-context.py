"""Check explicitly new junction faces against every active bounded native neighbor."""
from pathlib import Path
import argparse,hashlib,json,time,sys
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology,winding_contains
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/maxillary-native-v1/junction-repair-v1';C=B/'candidate-v4'
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,default=lambda v:v.tolist()if hasattr(v,'tolist')else v.item())+'\n')
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--start',type=int,default=0);ap.add_argument('--count',type=int,default=30);a=ap.parse_args();start=time.monotonic();plan=json.loads((B/'plan.json').read_text());cp=ROOT/plan['source_catalogue'];assert sha(cp)==plan['source_master_sha256'];cat=json.loads(cp.read_text());hidden=set(cat.get('defaultHidden',[]));m=Mesh.load(C/'cohort.npz');an={k:v for k,v in np.load(C/'ancestry.npz').items()};newids=np.flatnonzero(an['source_face_ids']<0);patch=Mesh(m.vertices,m.faces[newids],m.frame,m.attributes,m.source);tt=patch.triangles;lo=tt.min(1);hi=tt.max(1);selected={r['id']for r in plan['sources']if r['source']in plan['selected_whole_source_indices']};active=[];near=[]
 for r in cat['parts']:
  if r['id']in hidden or r.get('suppressed')or r.get('defaultHidden')or r['id']in selected:continue
  assert 'bounds'in r;active.append(r['id']);b=np.array(r['bounds']);hits=np.all(b[0]<=hi+2e-6,axis=1)&np.all(b[1]>=lo-2e-6,axis=1)
  if np.any(hits):near.append(r)
 O=C/f'independent-fixed-context-{a.start:03}-{min(a.start+a.count,len(near)):03}.json';rows=[];reader=AtlasReader();tree=BVH(patch)
 for index,r in list(enumerate(near))[a.start:a.start+a.count]:
  sm=reader.read_part(cat,r,frame=m.frame,catalogue_path=cp,catalogue_sha=plan['source_master_sha256']);sb=sm.bounds;mb=np.array(r['bounds']);bounds_valid=bool(np.all(sb[0]>=mb[0]-2e-6)and np.all(sb[1]<=mb[1]+2e-6))
  other=BVH(sm);p,s=pair_contacts(patch,sm,0.,first_bvh=tree,second_bvh=other,max_pairs=100000,max_tests=500000);assert s['complete_within_pair_budget'];wrong=np.isin(p['relation'],['crossing','coplanar-overlap']);pp=p['pairs'][wrong];t=topology(sm);inside=[];probe_count=0
  if t['closed_oriented_manifold']:
   pts=np.unique(np.concatenate([tt.reshape(-1,3),tt.mean(1)]),axis=0);pts=pts[np.all(pts>=sb[0],axis=1)&np.all(pts<=sb[1],axis=1)];probe_count=len(pts)
   if len(pts):
    w=winding_contains(sm,pts,boundary_tolerance_m=1e-10);inside=pts[np.array(w['inside'],bool)].tolist()
  rows.append(dict(index=index,id=r['id'],system=r.get('system'),geometry_sha256=sm.geometry_sha,attributes_sha256=sm.attributes_sha,actual_bounds_m=sb,metadata_bounds_contain_actual=bounds_valid,strict_crossings=len(pp),strict_witnesses=[dict(candidate_face=int(newids[x]),fixed_face=int(y),candidate_triangle_m=patch.triangles[x].tolist(),fixed_triangle_m=sm.triangles[y].tolist())for x,y in pp[:12]],contact_summary=s,closed_oriented_manifold=t['closed_oriented_manifold'],inside_probe_count=len(inside),tested_containment_probes=probe_count,inside_probe_examples_m=inside[:12]))
  report=dict(schema='maxillary-v4-active-fixed-context/v1',source_master_sha256=plan['source_master_sha256'],candidate_sha256=sha(C/'cohort.npz'),new_junction_faces=len(newids),active_fixed_records=len(active),metadata_broad_phase_neighbors=len(near),metadata_no_bounds=0,start=a.start,count=len(rows),expected_stage_count=min(a.count,max(0,len(near)-a.start)),rows=rows,complete=len(rows)==min(a.count,max(0,len(near)-a.start)),elapsed_seconds=time.monotonic()-start,scope='All active master metadata bounds examined, actual candidate new faces tested against actual decoded nearby fixed triangles. Original source-owned portions have unchanged positions and remain qualified source contacts. No peer exemptions or anatomy changes. Open tissue sheets do not permit volume containment classification.')
  report['passed']=report['complete']and not any(x['strict_crossings']or x['inside_probe_count']for x in rows);write(O,report)
 print(json.dumps(dict(file=str(O),complete=report['complete'],passed=report['passed'],neighbors=len(near),count=len(rows),strict=[(x['id'],x['strict_crossings'],x['inside_probe_count'])for x in rows if x['strict_crossings']or x['inside_probe_count']],elapsed=report['elapsed_seconds'])))
if __name__=='__main__':main()
