"""Exact new-inner witnesses and source stations; no fitting."""
from pathlib import Path
import json,sys,time
from collections import defaultdict
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';C=P/'three-course-finite-educational-v1';I=P/'temporal-source-incidence-v2';OUT=I/'inner-crossing-witnesses-v1';FRAME='model:local-ta98:assembled-metres';NEW={127,142,143,248,249,295}

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(F/'three-course-source-charts.npz');t=np.load(I/'complete-source-charts.npz');q=np.load(C/'assembly-v1/wall-with-source-and-derived-corner-normals.npz');generation=json.loads((C/'generation-complete.json').read_text());rows={r['original_id']:r for r in generation['classifiers']};allowed={'M3':set(z['M3_root_patch_face_ids'].tolist()),'main':set(z['main_root_to_M3_patch_face_ids'].tolist())};mask=[];ancestry=[];labels={}
 for name,vkey,rkey in [('M3','M3_vertices','M3_rings'),('main','branch_vertices','main_rings')]:
  lookup=np.full(len(z[vkey]),-1,int)
  for station,ring in enumerate(z[rkey]):lookup[ring]=station
  labels[name]=lookup
 for i,(oid,fid)in enumerate(zip(q['original_id'],q['source_face_id'])):
  row=rows[int(oid)];course=row['course']
  if course not in allowed or int(fid)>=row['source_faces']:continue
  sid=row['source_face_ids'][int(fid)]
  if sid not in allowed[course]:continue
  sf=z['M3_faces'][sid]if course=='M3'else z['branch_faces'][sid];mask.append(i);ancestry.append(dict(role=row['role'],course=course,original_wall_face=i,source_face=sid,source_stations=sorted(set(map(int,labels[course][sf])))))
 origin=np.array([0.,1.75,0.]);positions=(q['vertices']-origin).astype('<f4').astype(float)+origin;wall=Mesh(positions,q['faces'][mask],FRAME);temporal=Mesh.load(I/'complete-temporal-source.npz');arr,summary=pair_contacts(wall,temporal,0.);assert summary['complete_within_pair_budget'];pairs=arr['pairs'][np.isin(arr['relation'],['crossing','coplanar-overlap'])];temporal_labels=np.full(len(t['temporal_vertices']),-1,int)
 for station,ring in enumerate(t['component0_rings']):temporal_labels[ring]=station
 group=defaultdict(list)
 for first,second in pairs:
  a=ancestry[int(first)];r=dict(**a,local_wall_face=int(first),temporal_source_face=int(second),temporal_stations=sorted(set(map(int,temporal_labels[t['temporal_faces'][second]]))),new_peer_face=bool(int(second)in NEW),wall_triangle_m=wall.triangles[first].tolist(),temporal_triangle_m=temporal.triangles[second].tolist());group[(a['role'],bool(int(second)in NEW))].append(r)
 output=[]
 for (role,new),hits in group.items():output.append(dict(role=role,new_peer_face=new,pairs=len(hits),wall_stations=sorted(set(k for r in hits for k in r['source_stations'])),temporal_stations=sorted(set(k for r in hits for k in r['temporal_stations'])),temporal_source_faces=sorted(set(r['temporal_source_face']for r in hits))))
 atomic_json(OUT/'crossing-witnesses.json',dict(source_scope='Actual emitted local wall and complete temporal source part, all components retained.',summary=summary,groups=output,witnesses=[r for hits in group.values()for r in hits],new_peer_source_face_ids=sorted(NEW),geometry_modified=False,elapsed_seconds=time.monotonic()-start));atomic_npz(OUT/'actual-contact-pairs.npz',pairs=pairs,local_wall_to_validated_whole_face_ids=np.array(mask,int),frame=np.array(FRAME));print(json.dumps(dict(groups=output,elapsed_seconds=time.monotonic()-start)))

if __name__=='__main__':main()
