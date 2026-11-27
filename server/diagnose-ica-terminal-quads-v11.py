"""Read-only replay and deterministic diagonal quality diagnosis of fixed v10 geometry."""
from pathlib import Path
import sys,json,hashlib,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
B=ROOT/".local-models/geometry-workbench/circulatory-completion/anterior-cerebral-segments-v1";C=B/"native-carotid-terminal-cohort-v7b";O=B/"ica-mca-quality-triangulation-v11";sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,default=lambda v:v.tolist()if hasattr(v,"tolist")else v.item())+"\n")
def replay(side):
 plan=json.loads((C/"joint-interface-plan-v8.json").read_text());expanded=json.loads((B/"ica-mca-three-ring-expanded-plan-v10.json").read_text());source=Mesh.load(C/f"ica-{side}.npz");mca=Mesh.load(C/f"mca-{side}.npz");pp=next(x for x in plan["sides"]if x["side"]==side);ep=next(x for x in expanded["sides"]if x["side"]==side);v=source.vertices.copy();attrs={k:a.copy()for k,a in source.attributes.items()};endaxis=np.array(pp["mca_proximal_profile"]["best_plane_unit_normal"])
 for ids,angle in zip(ep["source_control_ids"][:2],ep["rotation_angles_deg"][:2]):
  ids=np.array(ids,int);centre=source.vertices[ids].mean(0);_,_,vh=np.linalg.svd(source.vertices[ids]-centre);axis=vh[-1]
  if np.dot(axis,endaxis)<0:axis=-axis
  x,y,z=axis;skew=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]]);t=np.deg2rad(angle);rotation=np.eye(3)+np.sin(t)*skew+(1-np.cos(t))*(skew@skew);v[ids]=(source.vertices[ids]-centre)@rotation.T+centre;attrs["normals"][ids]=np.rint(source.attributes["normals"][ids].astype(float)@rotation.T).astype(np.int16)
 ii=np.array(pp["correspondence"]["ica_distal_control_ids"],int);mi=np.array(pp["correspondence"]["mca_control_ids_in_ica_order"],int);v[ii]=mca.vertices[mi];attrs["normals"][ii]=mca.attributes["normals"][mi];q=Mesh(v,source.faces,source.frame,attrs)
 if side=="l":assert q.geometry_sha==json.loads((B/"ica-mca-seam-three-ring-v10/raw-l-surface-audit.json").read_text())["geometry_sha256"]
 rings,_=align_station_correspondence(source.vertices,source.faces,recover_station_rings(source.vertices,source.faces)[0]["rings"]);rings=np.array(rings)
 if np.linalg.norm(source.vertices[rings[0]].mean(0)-np.array(pp["ica_distal_profile"]["centre_m"]))<np.linalg.norm(source.vertices[rings[-1]].mean(0)-np.array(pp["ica_distal_profile"]["centre_m"])):rings=rings[::-1]
 return source,mca,q,rings

def cyclic_same(a,b):return any(np.array_equal(np.roll(a,k),b)for k in range(3))
def quality(v,ff,centres,station):
 p=v[ff];cross=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);norm=np.linalg.norm(cross,axis=1);radial=(v-centres[station])[ff].sum(1);radial/=np.linalg.norm(radial,axis=1)[:,None];signed=.5*np.einsum("ij,ij->i",cross,radial);align=2*signed/norm;area=.5*norm;edge2=np.square(p[:,1]-p[:,0]).sum(1)+np.square(p[:,2]-p[:,1]).sum(1)+np.square(p[:,0]-p[:,2]).sum(1);aspect=4*np.sqrt(3)*area/edge2
 return dict(projected_signed_areas_m2=signed,normal_radial_alignments=align,triangle_areas_m2=area,aspect_quality=aspect,valid_outward_nondegenerate=bool(np.all(signed>1e-18)and np.all(norm>2e-18)),score=[float(signed.min()),float(align.min()),float(aspect.min())])
def diagnose(side):
 source,mca,q,rings=replay(side);station=np.zeros(len(q.vertices),int);station[rings]=np.arange(len(rings))[:,None];centres=q.vertices[rings].mean(1);keys={tuple(sorted(t)):i for i,t in enumerate(source.faces)};rows=[]
 for interval in[7,8,9]:
  left,right=rings[interval],rings[interval+1]
  for k in range(12):
   j=(k+1)%12;a,b,c,d=map(int,[left[k],left[j],right[j],right[k]]);choices=[np.array([[a,c,b],[a,d,c]]),np.array([[a,d,b],[b,d,c]])];matches=[]
   for n,ff in enumerate(choices):
    if all(tuple(sorted(t))in keys for t in ff):matches.append((n,[keys[tuple(sorted(t))]for t in ff]))
   assert len(matches)==1;original,faceids=matches[0];orientation=1 if cyclic_same(choices[original][0],source.faces[faceids[0]])else-1
   if orientation==-1:choices=[ff[:,::-1]for ff in choices]
   qq=[quality(q.vertices,ff,centres,station)for ff in choices];valid=[i for i,x in enumerate(qq)if x["valid_outward_nondegenerate"]];selected=max(valid,key=lambda i:tuple(qq[i]["score"]))if valid else None;rows.append(dict(interval=interval,quad=k,corners=[a,b,c,d],corner_positions_m=q.vertices[[a,b,c,d]],source_face_ids=faceids,source_diagonal=original,source_faces=source.faces[faceids],candidate_diagonals=[dict(diagonal=i,faces=ff,quality=qc)for i,(ff,qc)in enumerate(zip(choices,qq))],selected_diagonal=selected,selected_by="Maximize minimum positive outward projected trianglearea, thenminimumradialnormalalignment, thenminimumaspectquality; no selfcounts used.",involves_failed_pair=any(i in[209,230]for i in faceids)))
 return q,rows
if __name__=="__main__":
 assert sha(ROOT/".local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json")=="1befb0177aa7f3209ece73dcf28bc03c05bd25cc429c6c5360982d2b62106114";O.mkdir(exist_ok=True);allrows=[]
 for side in["l","r"]:
  q,rows=diagnose(side);allrows.append(dict(side=side,v10_replayed_geometry_sha256=q.geometry_sha,v10_attributes_sha256=q.attributes_sha,quads=rows));print(side,"invalidquads",[(r["interval"],r["quad"])for r in rows if r["selected_diagonal"]is None],"changeddiagonals",sum(r["selected_diagonal"]is not None and r["selected_diagonal"]!=r["source_diagonal"]for r in rows),"failedpairquads",[(r["interval"],r["quad"],r["source_diagonal"],r["selected_diagonal"],[x["quality"]["score"]for x in r["candidate_diagonals"]])for r in rows if r["involves_failed_pair"]])
 write(O/"quad-diagnostic.json",dict(complete=True,geometry_changed=False,self_audit_run=False,sides=allrows))
