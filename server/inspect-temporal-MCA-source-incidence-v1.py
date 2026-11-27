"""Complete temporal-source component/port graph; no geometry modification."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.audit import topology
from geometry_tools.contact import pair_contacts
from geometry_tools.algorithms.tube_station_topology import recover_station_rings,align_station_correspondence
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';F=P/'secondary-port-feasibility-v1';C=P/'three-course-finite-educational-v1';CTX=C/'current-context-v1';OUT=P/'temporal-source-incidence-v2';FRAME='model:local-ta98:assembled-metres';PID='ZA:Temporal branches of middle cerebral artery.l';NEW=[127,142,143,248,249,295]

def chart(mesh):
 v,inv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inv[mesh.faces];groups=recover_station_rings(v,f);rows=[]
 for group in groups:
  rings,phase=align_station_correspondence(v,f[group['face_ids']],group['rings']);rings=np.array(rings);rows.append({**group,'rings':rings,'centers':v[rings].mean(1),'phase':phase})
 return v,f,rows

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);plan=json.loads((CTX/'scan-plan.json').read_text());source=CTX/'source-atlas.json';assert digest_file(source)==plan['source_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();mesh=reader.read_part(atlas,parts[PID],frame=FRAME);assert mesh.geometry_sha=='fd7f484f6a323ae30b9d4db26368512209aa3b6c27ab989412d2dca49a45b6a0';mesh.save(OUT/'complete-temporal-source.npz');v,f,groups=chart(mesh);z=np.load(F/'three-course-source-charts.npz');known=[('M3',z['M3_vertices'],z['M3_faces'],z['M3_rings']),('main-postcentral',z['branch_vertices'],z['branch_faces'][z['main_source_face_ids']],z['main_rings']),('secondary-postcentral',z['branch_vertices'],z['branch_faces'][z['twig_source_face_ids']],z['twig_rings'])];charts=[];candidateports=[];arrays=dict(temporal_vertices=v,temporal_faces=f,frame=np.array(FRAME));ownership=np.full(len(f),-1,int)
 for index,group in enumerate(groups):
  ids=group['face_ids'];ownership[ids]=index;rings=group['rings'];centers=group['centers'];endpoints=[]
  for station in [0,len(rings)-1]:
   ring=v[rings[station]];point=centers[station];radius=np.linalg.norm(ring-point,axis=1);matches=[]
   for name,pv,pf,pr in known:
    pc=pv[pr].mean(1);dist=np.linalg.norm(pc-point,axis=1);closest=int(dist.argmin());parentring=pv[pr[closest]];peraxis=.5*(np.mean(abs(np.spacing(ring.astype('<f4')).astype(float)),axis=0)+np.mean(abs(np.spacing(parentring.astype('<f4')).astype(float)),axis=0));roundoff=float(np.linalg.norm(peraxis)+64*np.finfo(float).eps*max(np.max(abs(ring)),np.max(abs(parentring))));segments=pc[1:]-pc[:-1];t=np.clip(np.einsum('ij,ij->i',point-pc[:-1],segments)/np.einsum('ij,ij->i',segments,segments),0,1);near=pc[:-1]+segments*t[:,None];sd=np.linalg.norm(near-point,axis=1);segment=int(sd.argmin());row=dict(course=name,nearest_station=closest,station_center_residual_m=float(dist[closest]),derived_two_ring_Float32_roundoff_bound_m=roundoff,station_match_within_source_storage_roundoff=bool(dist[closest]<=roundoff),nearest_centerline_segment=segment,segment_parameter=float(t[segment]),centerline_segment_distance_m=float(sd[segment]));matches.append(row)
    if row['station_match_within_source_storage_roundoff']:candidateports.append(dict(temporal_component=index,temporal_station=station,**row))
   endpoints.append(dict(station=station,center_m=point.tolist(),mean_outer_radius_m=float(radius.mean()),source_ring_vertex_ids=rings[station].tolist(),known_course_matches=matches))
  face_rows=[]
  for fid in NEW:
   if fid in set(ids.tolist()):face_rows.append(dict(source_face=fid,component=index,stations=sorted(set(map(int,group['station_of_vertex'][f[fid]])))))
  charts.append(dict(component=index,source_component=int(group['source_component']),source_faces=len(ids),source_face_min=int(ids.min()),source_face_max=int(ids.max()),stations=len(rings),ring_size=len(rings[0]),ports=endpoints,newly_crossed_source_faces=face_rows));arrays['component'+str(index)+'_rings']=rings;arrays['component'+str(index)+'_source_face_ids']=ids
 assert np.all(ownership>=0);atomic_npz(OUT/'complete-source-charts.npz',**arrays)
 # Whole-source contact attribution, not just the cropped six new witnesses.
 contacts=[]
 for name,pv,pf,pr in known:
  arr,summary=pair_contacts(mesh,Mesh(pv,pf,FRAME),0.);assert summary['complete_within_pair_budget'];pairs=arr['pairs'][np.isin(arr['relation'],['crossing','coplanar-overlap'])];bycomponent=[]
  for index,group in enumerate(groups):
   selected=pairs[ownership[pairs[:,0]]==index]
   if not len(selected):continue
   stations=group['station_of_vertex'][f[selected[:,0]]];bycomponent.append(dict(component=index,pairs=len(selected),temporal_station_range=[int(stations.min()),int(stations.max())],temporal_source_face_ids=sorted(set(map(int,selected[:,0]))),known_course_local_face_ids=sorted(set(map(int,selected[:,1])))))
  contacts.append(dict(known_course=name,summary=summary,crossing_components=bycomponent));atomic_npz(OUT/(name+'-contact-pairs.npz'),pairs=pairs)
 inventory=json.loads((B/'source-bound-inventory.json').read_text());authority=[dict(id=r['authority']['id'],name_en=r['authority']['name_en'],name_la=r['authority']['name_la'],parent_id=r['authority']['parent_id'])for r in inventory['terms']];qualifier='Native temporal branches source collection only; no individual TA98 cortical leaf assignment or complete .061 collective coverage is inferred from a geometric component.';report=dict(schema='temporal-MCA-native-source-incidence/v1',installed=False,geometry_modified=False,current_source_sha256=plan['source_sha256'],part_id=PID,geometry_sha256=mesh.geometry_sha,attributes_sha256=mesh.attributes_sha,part_ta98Term=parts[PID].get('ta98Term'),concept_id=parts[PID]['conceptId'],concept_ta98Term=next(c for c in atlas['concepts']if c['id']==parts[PID]['conceptId']).get('ta98Term'),full_source_faces=len(f),source_components=charts,newly_crossed_faces=NEW,newly_crossed_faces_covered_exactly_once=all(np.sum([fid in set(g['face_ids'].tolist())for g in groups])==1 for fid in NEW),endpoint_station_matches=candidateports,whole_known_course_contacts=contacts,authority_inventory=authority,term_qualification=qualifier,source_cortex_crossings_remain_separate_unaccepted=True,construction_authorized_by_this_report=False,elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'incidence-proof.json',report);print(json.dumps(dict(full_source_faces=len(f),components=[dict(component=r['component'],faces=r['source_faces'],stations=r['stations'],new_face_stations=r['newly_crossed_source_faces'])for r in charts],endpoint_station_matches=candidateports,contact_counts=[dict(course=r['known_course'],crossing_pairs=r['summary']['crossing_face_pairs'])for r in contacts],elapsed_seconds=report['elapsed_seconds'])))

if __name__=='__main__':main()
