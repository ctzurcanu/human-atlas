"""Read-only localization of source c1/infra-orbital wall-origin overlap."""
from pathlib import Path
from collections import Counter
import sys,json,hashlib,numpy as np,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology,winding_contains
from geometry_tools.distance import barycentric
B=ROOT/".local-models/geometry-workbench/circulatory-completion/maxillary-native-v1/junction-repair-v1";P=json.loads((B/"plan.json").read_text());COURSES=np.load(B/"source-course-arrays.npz");sources={x["source"]:x for x in P["sources"]};sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def close(m):
 f=m.faces;ee=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);count=Counter(tuple(sorted(x))for x in ee);unused={tuple(e)for e in ee if count[tuple(sorted(e))]==1};v=list(m.vertices);ff=[]
 while unused:
  edge=next(iter(unused));loop=[edge[0]];here=edge[0]
  while True:
   e=next(e for e in unused if e[0]==here);unused.remove(e);here=e[1]
   if here==loop[0]:break
   loop.append(here)
  c=len(v);v.append(m.vertices[loop].mean(0));ff.extend([[loop[(j+1)%len(loop)],loop[j],c]for j in range(len(loop))])
 q=Mesh(np.array(v),np.r_[f,np.array(ff)],m.frame);assert topology(q)["closed_oriented_manifold"];return q

def run():
 assert sha(ROOT/P["source_catalogue"])==P["source_master_sha256"];started=time.monotonic();out=[]
 for sid,pid,side in[(0,10,"l"),(1,11,"r")]:
  s=Mesh.load(sources[sid]["path"]);p=Mesh.load(sources[pid]["path"]);r=COURSES[f"s{sid}_c1_rings"];pr=COURSES[f"s{pid}_c0_rings"];station=np.full(len(s.vertices),-1,int);station[r]=np.arange(len(r))[:,None];pstation=np.full(len(p.vertices),-1,int);pstation[pr]=np.arange(len(pr))[:,None];fids=np.flatnonzero(station[s.faces[:,0]]>=0);child=Mesh(s.vertices,s.faces[fids],s.frame,s.attributes);cp,summary=pair_contacts(child,p,0.,max_pairs=250000,max_tests=1000000);bad=np.isin(cp["relation"],["crossing","coplanar-overlap"]);pairs=cp["pairs"][bad];sourcepairs=np.c_[fids[pairs[:,0]],pairs[:,1]];cs=station[s.faces[sourcepairs[:,0]]];ps=pstation[p.faces[sourcepairs[:,1]]];c=s.vertices[r].mean(1);pc=p.vertices[pr].mean(1);d=np.diff(pc,axis=0);t=np.clip(np.einsum("ij,ij->i",c[0]-pc[:-1],d)/np.einsum("ij,ij->i",d,d),0,1);q=pc[:-1]+t[:,None]*d;dd=np.linalg.norm(q-c[0],axis=1);j=int(dd.argmin());radial=(c[0]-q[j])/dd[j];tangent=c[1]-c[0];tangent/=np.linalg.norm(tangent);classifier=close(p);queries=np.r_[s.vertices[r[:9]].reshape(-1,3),c[:9]];w=winding_contains(classifier,queries);states=[];ptree=BVH(p)
  for i in range(9):
   ds=np.linalg.norm(pc[:-1]+np.clip(np.einsum("ij,ij->i",c[i]-pc[:-1],d)/np.einsum("ij,ij->i",d,d),0,1)[:,None]*d-c[i],axis=1);closest=ptree.closest_point(c[i]);states.append(dict(station=i,centre_m=c[i].tolist(),centre_axis_distance_m=float(ds.min()),centre_surface_distance_m=closest["distance_m"],profile_controls_inside_parent=sum(w["inside"][i*12:(i+1)*12]),centre_inside_parent=w["inside"][108+i],ambiguous_profile_controls=sum(w["boundary_ambiguous"][i*12:(i+1)*12])))
  np.savez_compressed(B/f"{side}-c1-infra-original-crossings.npz",source_face_pairs=sourcepairs,child_station_per_corner=cs,parent_station_per_corner=ps,child_centres=c,parent_centres=pc);out.append(dict(side=side,child_source=sources[sid]["id"],child_component=1,parent_source=sources[pid]["id"],source_crossing_pairs=int(bad.sum()),child_crossing_strips=sorted(set((int(x.min()),int(x.max()))for x in cs)),parent_crossing_strips=sorted(set((int(x.min()),int(x.max()))for x in ps)),child_crossing_station_min=int(cs.min()),child_crossing_station_max=int(cs.max()),parent_crossing_station_min=int(ps.min()),parent_crossing_station_max=int(ps.max()),root_parent_axis_segment=j,root_parent_axis_fraction=float(t[j]),root_axis_gap_m=float(dd[j]),initial_tangent_dot_away_from_parent_axis=float(tangent@radial),root_profile_mean_radius_m=float(np.linalg.norm(s.vertices[r[0]]-c[0],axis=1).mean()),parent_interpolated_mean_radius_m=float((1-t[j])*np.linalg.norm(p.vertices[pr[j]]-pc[j],axis=1).mean()+t[j]*np.linalg.norm(p.vertices[pr[j+1]]-pc[j+1],axis=1).mean()),profiles=states,scope="Geometrysupportsorrefuteslocalwall-origintakeoffonly;noclinicalbranchidentityorcodeassigned."))
 report=dict(complete=True,source_master_sha256=P["source_master_sha256"],source_plan_sha256=sha(B/"plan.json"),geometry_modified=False,sides=out,elapsed_seconds=time.monotonic()-started);(B/"c1-infra-root-diagnostic-v1.json").write_text(json.dumps(report,indent=2));print(json.dumps(report))
if __name__=="__main__":run()
