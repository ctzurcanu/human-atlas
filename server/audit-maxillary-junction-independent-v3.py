"""Independent literal-bank/source-ancestry maxillary exterior gate. No fitting or mutation."""
from pathlib import Path
from collections import Counter,defaultdict
import sys,json,hashlib,time,numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.model import Mesh
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from geometry_tools.audit import topology
from geometry_tools.distance import point_triangle
B=ROOT/".local-models/geometry-workbench/circulatory-completion/maxillary-native-v1/junction-repair-v1";C=B/"candidate-v3";sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()
def write(p,x):Path(p).write_text(json.dumps(x,indent=2,default=lambda v:v.tolist()if hasattr(v,"tolist")else v.item())+"\n")
def run():
 start=time.monotonic();plan=json.loads((B/"plan.json").read_text());proof=json.loads((C/"proof.json").read_text());assert proof["complete"]and sha(ROOT/plan["source_catalogue"] )==plan["source_master_sha256"]==proof["source_master_sha256"];assert sha(C/"cohort.npz")==proof["candidate_sha256"];m=Mesh.load(C/"cohort.npz");an=np.load(C/"ancestry.npz");course=np.load(B/"source-course-arrays.npz");sources={int(x["source"]):x for x in plan["sources"]};component_maps={};source_meshes={};retained=[x for x in plan["fixed_terminal_ports"]if not x["root_matches"]];consumed=[x for x in plan["fixed_terminal_ports"]if x["root_matches"]];assert len(retained)==24 and len(consumed)==20;expected_edges=set();retained_keys=[]
 for sid in plan["selected_whole_source_indices"]:
  row=sources[sid];assert sha(row["path"])==row["sha256"];sm=Mesh.load(row["path"]);source_meshes[sid]=sm;cm=np.full(len(sm.faces),-1,int)
  for cr in row["courses"]:
   rings=course[cr["prefix"]+"_rings"];owned=np.isin(sm.faces[:,0],rings.ravel());cm[owned]=cr["component"]
  assert np.all(cm>=0);component_maps[sid]=cm
 for port in retained:
  sm=source_meshes[port["source"]];rings=course[port["ring_array"]];ids=rings[port["ring_station"]];v=sm.vertices[ids];retained_keys.append(dict(source=port["source"],component=port["component"],end=port["end"],controls=len(ids)))
  expected_edges|={frozenset([tuple(v[j]),tuple(v[(j+1)%len(v)])])for j in range(len(v))}
 verts,inv=np.unique(m.vertices,axis=0,return_inverse=True);f=inv[m.faces];edges=np.concatenate([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);faceowners=np.tile(np.arange(len(f)),3);unique,first,counts=np.unique(np.sort(edges,axis=1),axis=0,return_index=True,return_counts=True);boundary=unique[counts==1];owners=faceowners[first[counts==1]];actual={frozenset(map(tuple,verts[e]))for e in boundary};missing=expected_edges-actual;extra=actual-expected_edges;exids=[j for j,e in enumerate(boundary)if frozenset(map(tuple,verts[e]))in extra];byrole=Counter();examples=[]
 for j in exids:
  fid=int(owners[j]);sid=int(an["source_indices"][fid]);sfid=int(an["source_face_ids"][fid]);comp=int(component_maps[sid][sfid])if sfid>=0 else -1;byrole[(sid,comp)]+=1
  if len(examples)<30:examples.append(dict(candidate_face_id=fid,source=sid,component=comp,source_face_id=sfid,edge_coordinates_m=verts[boundary[j]].tolist()))
 captri=[];capowners=[]
 for ex in proof["root_extensions"]:
  sid=ex["source"];comp=ex["component"];end=ex["end"];rr=course[f"s{sid}_c{comp}_rings"];ring=source_meshes[sid].vertices[rr[0 if end==0 else -1]]+np.array(ex["translation_m"]);centre=ring.mean(0)
  for k in range(len(ring)):captri.append([centre,ring[k],ring[(k+1)%len(ring)]]);capowners.append((sid,comp,end))
 captri=np.array(captri);capcounts=Counter();caperrors=[];unmapped=[]
 for j in exids:
  ep=verts[boundary[j]];d0,_,_=point_triangle(np.broadcast_to(ep[0],(len(captri),3)),captri);d1,_,_=point_triangle(np.broadcast_to(ep[1],(len(captri),3)),captri);distance=np.sqrt(np.maximum(d0,d1));q=int(distance.argmin());caperrors.append(float(distance[q]))
  if distance[q]<=1e-10:capcounts[capowners[q]]+=1
  else:unmapped.append(dict(edge=ep.tolist(),nearest_cap=capowners[q],maximum_endpoint_distance_to_same_cap_facet_m=float(distance[q])))
 expected_prefix=np.concatenate([source_meshes[sid].vertices for sid in plan["selected_whole_source_indices"]]);assert np.array_equal(an["source_position_prefix"],expected_prefix)and np.array_equal(m.vertices[:len(expected_prefix)],expected_prefix)
 for field in["normals","uvs"]:assert np.array_equal(m.attributes[field][:len(expected_prefix)],np.concatenate([source_meshes[sid].attributes[field]for sid in plan["selected_whole_source_indices"]]))
 tree=BVH(m);pairs,summary=pair_contacts(m,m,0.,first_bvh=tree,second_bvh=tree,selfcheck=True,max_pairs=600000,max_tests=2000000);assert summary["complete_within_pair_budget"];bad=np.isin(pairs["relation"],["crossing","coplanar-overlap"]);pp=pairs["pairs"][bad];by_source=Counter();by_component=Counter();witness=[]
 for pair in pp:
  desc=[]
  for fid in pair:
   sid=int(an["source_indices"][fid]);sfid=int(an["source_face_ids"][fid]);desc.append((sid,int(component_maps[sid][sfid])if sfid>=0 else -1,sfid))
  by_source[tuple(sorted([d[0]for d in desc]))]+=1;by_component[tuple(sorted([(d[0],d[1])for d in desc]))]+=1
  if len(witness)<40:witness.append(dict(candidate_faces=pair.tolist(),original_source_component_face=[list(d)for d in desc]))
 np.savez_compressed(C/"independent-boundary-self-witnesses.npz",geometric_boundary_edges=boundary,geometric_boundary_vertices=verts,boundary_adjacent_faces=owners,unexpected_boundary_indices=np.array(exids),strict_self_pairs=pp)
 report=dict(schema="maxillary-exterior-independent-failure-gate/v1",complete=True,passed=False,source_master_sha256=plan["source_master_sha256"],plan_sha256=sha(B/"plan.json"),candidate_sha256=sha(C/"cohort.npz"),candidate_triangles=len(m.faces),topology=topology(m),expected_external_ports=24,expected_external_edges=len(expected_edges),actual_geometric_boundary_edges=len(actual),missing_external_edges=len(missing),unexpected_geometric_boundary_edges=len(extra),unexpected_boundary_by_source_component=[dict(source=k[0],component=k[1],edges=v)for k,v in byrole.items()],unexpected_boundary_examples=examples,cap_derived_boundary_edges=[dict(source=k[0],component=k[1],end=k[2],edges=v)for k,v in capcounts.items()],unmapped_unexpected_edge_count=len(unmapped),unmapped_unexpected_examples=unmapped[:20],maximum_edge_to_nearest_cap_facet_distance_m=max(caperrors),original_source_position_normal_UV_prefix_literal=True,self_summary=summary,strict_crossings_by_original_source_pair=[dict(first=k[0],second=k[1],pairs=v)for k,v in by_source.items()],strict_crossings_by_original_component_pair=[dict(first=list(k[0]),second=list(k[1]),pairs=v)for k,v in by_component.items()],self_ancestry_examples=witness,retained_port_keys=retained_keys,elapsed_seconds=time.monotonic()-start,installable=False,scope="Actual Float64 source identity, exact-coordinate boundary and strict triangle-crossing audit; producerpassisnottrusted. No anatomy orcatalogue modified.");report["passed"]=not(missing or extra or len(pp));write(C/"independent-failure-audit.json",report);print(json.dumps({k:v for k,v in report.items()if k in["complete","passed","expected_external_edges","actual_geometric_boundary_edges","missing_external_edges","unexpected_geometric_boundary_edges","unexpected_boundary_by_source_component","cap_derived_boundary_edges","unmapped_unexpected_edge_count","strict_crossings_by_original_source_pair","strict_crossings_by_original_component_pair","elapsed_seconds"]}))
if __name__=="__main__":run()
