"""Endpoint-only incidence graph and explicit non-endpoint crossing witnesses."""
from pathlib import Path
import json,sys,time
import numpy as np
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh
from geometry_tools.contact import pair_contacts
B=ROOT/'.local-models/geometry-workbench/circulatory-completion/native-mca-controls-v3';P=B/'m3-postcentral-port-feasibility-v1';I=P/'temporal-source-incidence-v2';F=P/'secondary-port-feasibility-v1';OUT=I/'complete-graph-v1';FRAME='model:local-ta98:assembled-metres'

def nearest(point,port_ring,pv,pr):
 centers=pv[pr].mean(1);distance=np.linalg.norm(centers-point,axis=1);station=int(distance.argmin());parentring=pv[pr[station]];step=.5*(np.mean(abs(np.spacing(port_ring.astype('<f4')).astype(float)),axis=0)+np.mean(abs(np.spacing(parentring.astype('<f4')).astype(float)),axis=0));bound=float(np.linalg.norm(step)+64*np.finfo(float).eps*max(np.max(abs(port_ring)),np.max(abs(parentring))));return dict(nearest_station=station,station_center_residual_m=float(distance[station]),source_storage_roundoff_bound_m=bound,station_match=bool(distance[station]<=bound))

def main():
 start=time.monotonic();OUT.mkdir(exist_ok=False);z=np.load(I/'complete-source-charts.npz');known=np.load(F/'three-course-source-charts.npz');v=z['temporal_vertices'];f=z['temporal_faces'];courses=[('M3',known['M3_vertices'],known['M3_rings']),('main-postcentral',known['branch_vertices'],known['main_rings']),('secondary-postcentral',known['branch_vertices'],known['twig_rings'])]+[('temporal'+str(i),v,z['component'+str(i)+'_rings'])for i in [0,1]];edges=[];ports=[]
 for name,pv,pr in courses:
  for station in [0,len(pr)-1]:
   ring=pv[pr[station]];point=ring.mean(0);matches=[]
   for target,tv,tr in courses:
    if name==target:continue
    row=dict(target=target,**nearest(point,ring,tv,tr));matches.append(row)
    if row['station_match']:edges.append(dict(source=name,endpoint_station=station,source_endpoint_center_m=point.tolist(),**row))
   ports.append(dict(course=name,station=station,center_m=point.tolist(),mean_outer_radius_m=float(np.linalg.norm(ring-point,axis=1).mean()),comparisons=matches))
 # Exact existing temporal component-to-component surface contacts and their
 # source station ranges distinguish the genuine short-course root.
 a=z['component0_source_face_ids'];b=z['component1_source_face_ids'];arr,summary=pair_contacts(Mesh(v,f[a],FRAME),Mesh(v,f[b],FRAME),0.);assert summary['complete_within_pair_budget'];pairs=arr['pairs'][np.isin(arr['relation'],['crossing','coplanar-overlap'])];labels=[]
 for i in [0,1]:
  station=np.full(len(v),-1,int)
  for j,ring in enumerate(z['component'+str(i)+'_rings']):station[ring]=j
  labels.append(station)
 internal=dict(summary=summary,temporal0_station_range=[int(labels[0][f[a[pairs[:,0]]]].min()),int(labels[0][f[a[pairs[:,0]]]].max())]if len(pairs)else None,temporal1_station_range=[int(labels[1][f[b[pairs[:,1]]]].min()),int(labels[1][f[b[pairs[:,1]]]].max())]if len(pairs)else None);atomic_npz(OUT/'temporal-component-contact-pairs.npz',pairs=np.column_stack([a[pairs[:,0]],b[pairs[:,1]]]))
 incidence=json.loads((I/'incidence-proof.json').read_text());overlap=[]
 for name,pv,pr in courses[:3]:
  cp=np.load(I/(name+'-contact-pairs.npz'))['pairs'];pf=known['M3_faces']if name=='M3'else known['branch_faces'][known['main_source_face_ids']]if name=='main-postcentral'else known['branch_faces'][known['twig_source_face_ids']];station=np.full(len(pv),-1,int)
  for j,ring in enumerate(pr):station[ring]=j
  for i in [0,1]:
   own=set(z['component'+str(i)+'_source_face_ids'].tolist());selected=cp[np.array([int(k)in own for k in cp[:,0]],dtype=bool)]
   if not len(selected):continue
   trange=[int(labels[i][f[selected[:,0]]].min()),int(labels[i][f[selected[:,0]]].max())];krange=[int(station[pf[selected[:,1]]].min()),int(station[pf[selected[:,1]]].max())];overlap.append(dict(temporal_component=i,other_course=name,crossing_pairs=len(selected),temporal_station_range=trange,other_station_range=krange,contains_temporal_endpoint=trange[0]==0 or trange[1]==len(z['component'+str(i)+'_rings'])-1))
 witnesses=[]
 for fid in incidence['newly_crossed_faces']:
  stations=sorted(set(map(int,labels[0][f[fid]])));detail=[]
  for station in stations:
   ring=v[z['component0_rings'][station]];point=ring.mean(0);detail.append(dict(temporal_station=station,center_m=point.tolist(),known_course_matches=[dict(target=name,**nearest(point,ring,pv,pr))for name,pv,pr in courses[:3]]))
  witnesses.append(dict(source_face=fid,temporal_component=0,source_station_range=stations,at_any_temporal_endpoint=False,stations=detail))
 authority=incidence['authority_inventory'];contract=dict(schema='source-incidence-expansion-contract-only/v1',construction_authorized=False,installed=False,whole_courses_required=[dict(name=name,stations=len(pr),source_part_id='ZA:Temporal branches of middle cerebral artery.l'if name.startswith('temporal')else'ZA:Postcentral arterial branch.l'if 'postcentral'in name else'ZA:Middle cerebral artery (M3 segment).l')for name,_,pr in courses],declared_endpoint_graph=edges,undeclared_through_crossings=[r for r in overlap if r['other_course']=='main-postcentral'],scope='All complete source courses retained, no component assigned a falsely precise TA98 leaf. Endpoint graph alone does not authorize union across unrelated body crossings.',required_before_construction=['Root review of actual source incidence and rendered clash witnesses.','Separate genuine endpoint attachment from the temporal0/main-postcentral through-course crossings; blindly unioning all overlaps would invent arterial communication.','Keep both temporal source components, all fixed brain/skeleton/context, and original source calibre. Any later course correction needs a separately declared local domain and recipient check.','Native cortical source intersections remain separately unaccepted.'],term_qualification=dict(source_part_concept_term=incidence['concept_ta98Term'],source_part_term=incidence['part_ta98Term'],authority_names=authority,individual_temporal_leaf_identity_not_established=True));atomic_json(OUT/'expansion-contract.json',contract);report=dict(schema='complete-temporal-source-endpoint-graph/v1',current_source_sha256=incidence['current_source_sha256'],source_geometry_sha256=incidence['geometry_sha256'],all_source_components_and_ports_retained=True,endpoint_graph=edges,all_ports=ports,internal_temporal_contacts=internal,whole_crossing_station_ranges=overlap,newly_crossed_faces=witnesses,source_collection_term=incidence['concept_ta98Term'],interpretation='The actual temporal0 endpoint attaches at M3 station24. The six new wall-crossed faces are on temporal0 stations5–13, away from any endpoint; no exact endpoint incidence authorizes joining them to the existing M3-station36/postcentral root.',no_geometry_modified=True,no_construction_or_installation=True,elapsed_seconds=time.monotonic()-start);atomic_json(OUT/'graph-proof.json',report);print(json.dumps(dict(endpoint_graph=edges,internal_temporal_contacts={k:internal[k]for k in ['temporal0_station_range','temporal1_station_range']},whole_crossing_station_ranges=overlap,elapsed_seconds=report['elapsed_seconds'])))

if __name__=='__main__':main()
