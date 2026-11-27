"""Read-only actual source carotid fork/sinus geometry diagnostic."""
from pathlib import Path
from collections import Counter
import sys,json,hashlib,time,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh,AtlasReader
from geometry_tools.audit import audit,topology
from geometry_tools.spatial import BVH,box_distance_squared
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
CAT=ROOT/".local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json";B=ROOT/".local-models/geometry-workbench/circulatory-completion/carotid-bifurcation-sinus-v1c";EXPECTED="75aa4b959fdfa1500909b8b86cea8bdfc23146be70a78125b62310cd01c8c977";sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,default=lambda v:v.tolist()if hasattr(v,"tolist")else v.item())+"\n")
def ports(m):
 used=np.unique(m.faces);vertices,first,inverse=np.unique(m.vertices[used],axis=0,return_index=True,return_inverse=True);representative=used[first];mapping=np.full(len(m.vertices),-1,int);mapping[used]=inverse;f=mapping[m.faces];ee=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);counts=Counter(tuple(sorted(e))for e in ee);unused={tuple(e)for e in ee if counts[tuple(sorted(e))]==1};out=[]
 while unused:
  edge=next(iter(unused));here=edge[0];loop=[here]
  while True:
   choices=[e for e in unused if e[0]==here]
   if len(choices)!=1:raise ValueError("Non-simple directed source boundary")
   e=choices[0];unused.remove(e);here=e[1]
   if here==loop[0]:break
   loop.append(here)
  loop=representative[np.array(loop,int)].tolist();q=m.vertices[loop];c=q.mean(0);r=np.linalg.norm(q-c,axis=1);out.append(dict(control_ids=loop,controls=len(loop),centre_m=c.tolist(),radius_mean_m=float(r.mean()),radius_min_m=float(r.min()),radius_max_m=float(r.max())))
 return out

def run():
 assert sha(CAT)==EXPECTED;assert not B.exists();B.mkdir();started=time.monotonic();a=json.loads(CAT.read_text());reader=AtlasReader();frame="male-native-world-metres";meshes={};source_rows=[];sides=[];report=dict(complete=False,source_catalogue_sha256=EXPECTED,source_catalogue_path=str(CAT),sources=source_rows,sides=sides,geometry_mutated=False)
 for side,word in[("l","Left"),("r","Right")]:
  ids=dict(CCA=f"ZA:{word} common carotid artery",ECA=f"ZA:External carotid artery.{side}",ICA=f"TA98PART:ZA:Internal carotid artery.{side}:A12.2.06.002");srow=dict(side=side,source_ids=ids,interfaces=[]);sides.append(srow)
  for key,pid in ids.items():
   rec=next(p for p in a["parts"]if p["id"]==pid);assert not rec.get("suppressed");m=reader.read_part(a,rec,frame=frame,catalogue_path=CAT);m.save(B/f"{side}-{key}.npz");meshes[(side,key)]=m;au=audit(m,expected_components=1,require_closed=False,max_pairs=250000,max_tests=2000000);pp=ports(m);finite_annuli=[]
   if not pp and rec.get("pairedTissue"):
    for gi,g in enumerate(rec.get("tissueDrawGroups",[])):
     if g["count"]!=72:continue
     faceids=np.arange(g["start"]//3,(g["start"]+g["count"])//3);am=Mesh(m.vertices,m.faces[faceids],frame);loops=ports(am);assert len(loops)==2;loops.sort(key=lambda q:q["radius_mean_m"],reverse=True);outer,inner=loops;outer.update(finite_annular_tissue_port=True,draw_group_index=gi,source_face_ids=faceids.tolist(),inner_lumen_profile=inner,neighbor=g.get("neighbor"));pp.append(outer);finite_annuli.append(dict(draw_group=gi,outer=outer,inner=inner))
   row=dict(side=side,role=key,id=pid,name=rec["name"],geometry_sha256=m.geometry_sha,attributes_sha256=m.attributes_sha,vertices=len(m.vertices),triangles=len(m.faces),bounds_m=m.bounds.tolist(),ports=pp,audit=au,pairedTissue=bool(rec.get("pairedTissue")),finite_annular_mouths=finite_annuli);source_rows.append(row);write(B/"source-diagnostic.json",report)
  cca=next(x for x in source_rows if x["side"]==side and x["role"]=="CCA");upper=max(cca["ports"],key=lambda p:p["centre_m"][1]);srow["CCA_upper_port"]=upper
  for key in["ICA","ECA"]:
   row=next(x for x in source_rows if x["side"]==side and x["role"]==key);row["fork_nearest_port_index"]=int(np.argmin([np.linalg.norm(np.array(p["centre_m"])-upper["centre_m"])for p in row["ports"]]))if row["ports"]else None;m=meshes[(side,key)]
   if key=="ICA":
    rings,ph=align_station_correspondence(m.vertices,m.faces,recover_station_rings(m.vertices,m.faces)[0]["rings"]);rings=np.array(rings);c=m.vertices[rings].mean(1)
    if np.linalg.norm(c[-1]-upper["centre_m"])<np.linalg.norm(c[0]-upper["centre_m"]):rings=rings[::-1];c=c[::-1]
    arc=np.r_[0,np.cumsum(np.linalg.norm(np.diff(c,axis=0),axis=1))];rr=np.linalg.norm(m.vertices[rings]-c[:,None],axis=2);np.savez_compressed(B/f"{side}-ICA-stations.npz",rings=rings,centres=c,arc_m=arc,radii_m=rr);row["calibre_profiles"]=[dict(station=i,arc_m=float(ar),mean_radius_m=float(rad.mean()),min_radius_m=float(rad.min()),max_radius_m=float(rad.max()),centre_m=centre.tolist())for i,(ar,rad,centre)in enumerate(zip(arc,rr,c))];row["bulb_evidence"]=dict(first15mm_max_mean_radius_m=float(rr[arc<=.015].mean(1).max()),middle_third_median_mean_radius_m=float(np.median(rr[len(rr)//3:2*len(rr)//3].mean(1))),full_radius_min_m=float(rr.min()),full_radius_max_m=float(rr.max()),scope="Geometriccalibrevariationonly;doesnotestablishcarotidsinuswallhistologyorcorrectanatomicalbulbboundary.")
  for first,second in[("CCA","ECA"),("CCA","ICA"),("ICA","ECA")]:
   x=meshes[(side,first)];y=meshes[(side,second)];xb=BVH(x);yb=BVH(y);near=xb.closest_mesh(yb);cp,summary=pair_contacts(x,y,0.,first_bvh=xb,second_bvh=yb,max_pairs=250000,max_tests=2000000);np.savez_compressed(B/f"{side}-{first}-{second}-contacts.npz",**cp);xp=next(t for t in source_rows if t["side"]==side and t["role"]==first)["ports"];yp=next(t for t in source_rows if t["side"]==side and t["role"]==second)["ports"];pairs=[]
   for i,p in enumerate(xp):
    for j,q in enumerate(yp):
     v=x.vertices[p["control_ids"]];w=y.vertices[q["control_ids"]];pairs.append(dict(first_port=i,second_port=j,centre_gap_m=float(np.linalg.norm(np.array(p["centre_m"])-q["centre_m"])),literal_shared_controls=len(set(map(tuple,v))&set(map(tuple,w))),profile_control_Hausdorff_m=float(max(np.linalg.norm(v[:,None]-w[None],axis=2).min(1).max(),np.linalg.norm(w[:,None]-v[None],axis=2).min(1).max()))))
   pairs.sort(key=lambda q:q["centre_gap_m"]);srow["interfaces"].append(dict(first=first,second=second,whole_surface_nearest=near,contact=summary,closest_mouth_matches=pairs[:4]))
  rootpoint=np.array(upper["centre_m"]);lo=rootpoint-.035;hi=rootpoint+.035;context=[]
  for p in a["parts"]:
   if p.get("suppressed")or p.get("defaultHidden")or p["id"]in a["defaultHidden"]:continue
   structural=p.get("system")=="skeletal"or any(w in p["id"]for w in["CRICOID:male:","TA98HYOID:"])
   if not structural:continue
   bounds=np.array(p["bounds"])
   if box_distance_squared(lo,hi,bounds[0],bounds[1])>0:continue
   m=reader.read_part(a,p,frame=frame,catalogue_path=CAT);bt=BVH(m);q=bt.closest_point(rootpoint);context.append(dict(id=p["id"],name=p["name"],centre_to_surface_distance_m=q["distance_m"],nearest=q,geometry_sha256=m.geometry_sha,scope="CandidatefrommetadataROIthenactualtrianglequery;notglobalbodyskeletonminimum."))
  context.sort(key=lambda q:q["centre_to_surface_distance_m"]);srow["nearby_structural_context"]=context[:12];write(B/"source-diagnostic.json",report)
 assert sha(CAT)==EXPECTED;report.update(complete=True,elapsed_seconds=time.monotonic()-started,source_catalogue_unchanged=True);write(B/"source-diagnostic.json",report);print(json.dumps(dict(complete=True,elapsed_seconds=report["elapsed_seconds"],sources=[dict(id=x["id"],vertices=x["vertices"],faces=x["triangles"],components=len(x["audit"]["topology"]["components"]),mouths=len(x["ports"]),self=x["audit"]["self_intersections"]["defect_pairs"],calibre=x.get("bulb_evidence"))for x in source_rows],forks=[dict(side=x["side"],upper_port=x["CCA_upper_port"],interfaces=[dict(first=y["first"],second=y["second"],contact=y["contact"]["crossing_face_pairs"],mouthmatch=y["closest_mouth_matches"][0])for y in x["interfaces"]],context=[(y["id"],y["centre_to_surface_distance_m"])for y in x["nearby_structural_context"][:5]])for x in sides])))
if __name__=="__main__":run()
