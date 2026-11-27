"""Bounded fixed recipient-chart checks and literal port-normal interpolation."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.finite_tissue_boolean import barycentric
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';C=P/'finite-junction-educational-v1';S=C/'literal-splice-v1';OUT=C/'recipient-context-v1'

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);old=Mesh.load(S/'old-local-outer.npz');new=Mesh.load(S/'new-local-outer.npz');wall=Mesh.load(C/'intrinsic-v1/wall.npz');rows=[]
 # Only these actual frozen recipient charts are evaluated. No all-body
 # broad phase or proxy volume declares unrelated anatomy collision-free.
 for name in ['l-Parietal-lobe-1','l-Brain-Inner-2','l-Parietal-lobe-3','l-Interlobar-sulci-4']:
  peer=Mesh.load(B/(name+'.npz'));row=dict(recipient=name,fixed_geometry_sha256=peer.geometry_sha)
  for role,mesh in [('old_local_outer',old),('new_local_outer',new),('new_local_wall',wall)]:
   _,summary=pair_contacts(mesh,peer,0.);assert summary['complete_within_pair_budget'];row[role]=summary
  rows.append(row);atomic_json(OUT/(name+'-contacts.json'),row)
 q=np.load(S/'literal-splice-ancestry.npz');local=q['local_faces'];lineage=q['local_source_part_and_face'];corner=q['local_corner_interpolated_source_normal_words'];v=q['vertices'];charts=np.load(P/'source-port-charts.npz');sources=[Mesh.load(C/'parent-source-untouched.npz'),Mesh.load(C/'branch-source-untouched.npz')];normalrows=[]
 for role,index,key,station in [('parent',0,'parent_station_rings',33),('parent',0,'parent_station_rings',41),('branch',1,'branch_root_station_rings',176)]:
  points=charts[role+'_vertices'][charts[key][station]];portkeys={tuple(p)for p in points};tested=0;error=0.;geometric_error=0.;source=sources[index]
  for tri,ancestry,n in zip(local,lineage,corner):
   if ancestry[0]!=index:continue
   sourceface=source.faces[ancestry[1]];sq=source.vertices[sourceface]
   for p,normal in zip(v[tri],n):
    if tuple(p)not in portkeys:continue
    match=np.flatnonzero(np.all(sq==p,axis=1));assert len(match)==1;want=source.attributes['normals'][sourceface[match[0]]];error=max(error,float(np.max(abs(normal-want))));tested+=1
   # Every output side fragment retains its original planar face tangent.
   a=v[tri];b=sq;na=np.cross(a[1]-a[0],a[2]-a[0]);nb=np.cross(b[1]-b[0],b[2]-b[0]);na/=np.linalg.norm(na);nb/=np.linalg.norm(nb);geometric_error=max(geometric_error,float(np.linalg.norm(na-nb)))
  assert tested>=12 and error<1e-6 and geometric_error<1e-7
  normalrows.append(dict(role=role,station=station,port_corner_occurrences=tested,maximum_interpolated_source_normal_word_error=error,maximum_source_face_unit_normal_deviation=geometric_error,no_new_normal_or_tangent_discontinuity_at_literal_port=True))
 proof=dict(schema='local-M3-postcentral-recipient-context/v1',installed=False,recipients=rows,port_normals=normalrows,scope='Actual local patch versus four frozen nearby cerebral recipient charts only. No skeleton or context mesh moved. Secondary source-course crossings remain unresolved outside this local root; this does not establish vessel perfusion anatomy or whole-body packing.',elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'context-proof.json',proof);print(json.dumps(dict(recipients=[dict(name=r['recipient'],old=r['old_local_outer']['crossing_face_pairs'],new=r['new_local_outer']['crossing_face_pairs'],wall=r['new_local_wall']['crossing_face_pairs'])for r in rows],port_normals=normalrows,elapsed_seconds=proof['elapsed_seconds'])))

if __name__=='__main__':main()
