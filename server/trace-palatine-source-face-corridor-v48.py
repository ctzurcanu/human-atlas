"""Trace a literal SG face chain and measure its actual SG/SP medial interval."""
from pathlib import Path
import json,sys,time,heapq,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from collections import defaultdict
from geometry_tools.model import AtlasReader,Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';G=D/'global-palatine-actual-radius-v39';Q=D/'palatine-surface-portal-reachability-v45';O=D/'palatine-source-face-corridor-v48';assert not O.exists();O.mkdir();P=Q/'terminal-next-action-v46.json';gap=json.loads(P.read_bytes());reference=gap['four_target_geometric_gaps'][0];old=json.loads((G/'edge-search-evidence-v40.json').read_bytes());C=json.loads((G/'construction-contract.json').read_bytes());source=ROOT/C['source_catalogue'];assert digest_file(source)==C['source_catalogue_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';source_bind={p['part_id']:p for p in old['source_bindings']};sg=reader.read_part(atlas,parts['TA98TONGUE:male:styloglossus:left'],frame=frame,catalogue_path=source,catalogue_sha=C['source_catalogue_sha256']);sp=reader.read_part(atlas,parts['ZA:Stylopharyngeus muscle.l'],frame=frame,catalogue_path=source,catalogue_sha=C['source_catalogue_sha256']);assert sg.geometry_sha==source_bind['TA98TONGUE:male:styloglossus:left']['geometry_sha256'];assert sp.geometry_sha==source_bind['ZA:Stylopharyngeus muscle.l']['geometry_sha256'];a=np.array(reference['nearest_reachable_point_m']);b=np.array(reference['target_point_m']);_,face,near=igl.point_mesh_squared_distance(np.array([a,b]),sg.vertices,sg.faces);first,last=map(int,face);source_endpoints=near.copy();threshold=C['radius_m']+C['margin_m']
    contract=dict(schema='native-SG-SP-face-neighborhood-corridor-contract/v48',source_catalogue_sha256=C['source_catalogue_sha256'],terminal_gap_evidence_sha256=digest_file(P),source_count=179,source_profile_radius_m=C['radius_m'],clearance_margin_m=C['margin_m'],SG_geometry_sha256=sg.geometry_sha,SP_geometry_sha256=sp.geometry_sha,first_SG_face=first,last_SG_face=last,exact_endpoints_m=[a.tolist(),b.tolist()],
      method='Exact coordinate-welded SG triangle adjacency defines a continuous source-surface strip. Find its shortest face-centroid chain, cross literal shared-edge midpoints, subdivide those face-interior line segments at <=0.25 mm. For each source point project onto actual SP triangles and test the midpoint of that opposing pair plus fixed source endpoints against actual179-source context. Record SP correspondence discontinuities/clearance witnesses. This single face chain is a diagnostic, not proof that all SG/SP paths fail.',
      no_random_points_global_fields_or_radius_reduction=True,no_source_point_deformation=True,fixed_parent_station8=True,other3_routes_unchanged=C['unchanged_other_route_hashes'],internal_work_deadline_seconds=48,no_emitted_geometry_or_install=True)
    atomic_json(O/'construction-contract.json',contract)
    unique,inverse=np.unique(sg.vertices,axis=0,return_inverse=True);faces=inverse[sg.faces];owners=defaultdict(list)
    for i,f in enumerate(faces):
        for x,y in [(f[0],f[1]),(f[1],f[2]),(f[2],f[0])]:owners[tuple(sorted((int(x),int(y))))].append(i)
    adjacency=defaultdict(list);edge_map={};center=sg.triangles.mean(1);nonmanifold=[]
    for edge,ids in owners.items():
        if len(ids)!=2:
            if len(ids)>2:nonmanifold.append(dict(edge=list(edge),source_faces=ids))
            continue
        x,y=ids;mid=unique[list(edge)].mean(0);adjacency[x].append(y);adjacency[y].append(x);edge_map[tuple(sorted(ids))]=(edge,mid)
    distance={first:0.};previous={};queue=[(0.,first)];connected=False
    while queue:
        dd,f=heapq.heappop(queue)
        if dd!=distance[f]:continue
        if f==last:connected=True;break
        for t in adjacency[f]:
            nd=dd+float(np.linalg.norm(center[f]-center[t]))
            if nd<distance.get(t,np.inf):distance[t]=nd;previous[t]=f;heapq.heappush(queue,(nd,t))
    if not connected:atomic_json(O/'evidence.json',dict(completed=True,source_face_chain_exists=False,first_face=first,last_face=last,nonmanifold_edges=nonmanifold,no_diagnostic_guide=True,elapsed_seconds=time.monotonic()-began));return
    chain=[];f=last
    while f!=first:chain.append(f);f=previous[f]
    chain.append(first);chain=chain[::-1];waypoints=[source_endpoints[0]];transitions=[]
    for x,y in zip(chain[:-1],chain[1:]):
        edge,mid=edge_map[tuple(sorted((x,y)))];waypoints.append(mid);transitions.append(dict(first_face=x,second_face=y,shared_source_vertex_ids=list(edge),shared_edge_midpoint_m=mid.tolist()))
    waypoints.append(source_endpoints[1]);stations=[];source_faces=[]
    for i,(x,y)in enumerate(zip(waypoints[:-1],waypoints[1:])):
        count=max(1,int(np.ceil(np.linalg.norm(y-x)/.00025)))
        for k in range(count):stations.append(x+(y-x)*(k/count));source_faces.append(chain[i])
    stations.append(waypoints[-1]);source_faces.append(chain[-1]);stations=np.array(stations);d2,spface,spnear=igl.point_mesh_squared_distance(stations,sp.vertices,sp.faces);gaplen=np.sqrt(d2);medial=(stations+spnear)/2;proposed=np.vstack([a,medial,b]);sample=np.full(len(proposed),-1,dtype=int);sample[1:-1]=np.arange(len(stations));field=json.loads((G/'field-extension-evidence.json').read_bytes());query=np.load(G/'literal-surface-query-bank-v40.npz');assert digest_file(G/'literal-surface-query-bank-v40.npz')==old['query_bank_sha256'];tree=igl.AABB();tree.init(query['vertices'],query['faces']);squared,fi,closest=tree.squared_distance(query['vertices'],query['faces'],proposed);dist=np.sqrt(squared);safe=dist>=threshold;inside_witnesses=[];spec=importlib.util.spec_from_file_location('v48_source_occupancy',ROOT/'server/check-facial-cervical-feasibility.py');helper=importlib.util.module_from_spec(spec);spec.loader.exec_module(helper);complete=True;B=D.parent;plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes())
    for row in old['source_bindings']:
        if time.monotonic()-began>48:complete=False;break
        pid=row['part_id']
        if row['classifier']not in['literal_closed_source','temporary_native_conduit_end_caps_never_rendered']:continue
        if pid in parts:m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=C['source_catalogue_sha256']);allow=parts[pid].get('system')in{'arterial','venous'}
        elif pid=='CANDIDATE:glandular-0:passed-fixed-envelope':m=Mesh.load(ROOT/plan['extra_fixed_meshes'][0]['mesh']);allow=True
        else:
            name=pid.split(':')[1];z=np.load(D/(name+'-fields.npz'));m=Mesh(z['fields'][0],z['faces'],frame);allow=True
        assert m.geometry_sha==row['geometry_sha256'];cl,mode=helper.occupancy_classifier(m,allow);assert cl is not None;ids=np.flatnonzero(np.all(proposed>=cl.bounds[0],axis=1)&np.all(proposed<=cl.bounds[1],axis=1))
        if len(ids):inside=abs(igl.winding_number(cl.vertices,cl.faces,proposed[ids]))>.99;safe[ids[inside]]=False;inside_witnesses.extend(dict(point_index=int(i),part_id=pid,inside_actual_closed_source=True)for i in ids[inside])
        if pid in ['ZA:Facial artery.l','ZA:Palatopharyngeus muscle.l']:
            tt=igl.AABB();tt.init(m.vertices,m.faces);d2,_,_=tt.squared_distance(m.vertices,m.faces,proposed);dist=np.minimum(dist,np.sqrt(d2));safe&=np.sqrt(d2)>=threshold
    failures=[]
    for i in np.flatnonzero(~safe):
        faceid=int(fi[i]);owner=int(query['source_index'][faceid]);r=old['source_bindings'][owner];failures.append(dict(point_index=int(i),medial_point_m=proposed[i].tolist(),SG_station_index=int(sample[i]),SG_face=int(source_faces[sample[i]])if sample[i]>=0 else None,SP_face=int(spface[sample[i]])if sample[i]>=0 else None,opposing_surface_separation_m=float(gaplen[sample[i]])if sample[i]>=0 else None,limiting_source=r['part_id'],limiting_geometry_sha256=r['geometry_sha256'],source_face=int(query['native_face'][faceid]),surface_distance_m=float(np.sqrt(squared[i])),source_point_m=closest[i].tolist()))
    jumps=[]
    for i,(u,v)in enumerate(zip(spnear[:-1],spnear[1:])):
        if np.linalg.norm(v-u)>.001:jumps.append(dict(station_interval=i,SG_surface_step_m=float(np.linalg.norm(stations[i+1]-stations[i])),SP_projection_jump_m=float(np.linalg.norm(v-u)),first_SP_face=int(spface[i]),second_SP_face=int(spface[i+1])))
    atomic_npz(O/'actual-source-face-strip.npz',SG_surface_points=stations,SG_source_faces=np.array(source_faces),SP_nearest_points=spnear,SP_source_faces=spface,opposing_separation_m=gaplen,diagnostic_medial_points=proposed,point_surface_clearance_m=dist,point_pass=safe,frame=np.array(frame))
    result=dict(schema='native-SG-SP-continuous-face-strip-analysis/v48',source_catalogue_sha256=C['source_catalogue_sha256'],contract_sha256=digest_file(O/'construction-contract.json'),complete=complete,source_face_chain_exists=True,source_face_chain=chain,shared_edge_transitions=transitions,SG_source_surface_length_m=float(np.linalg.norm(np.diff(stations,axis=0),axis=1).sum()),source_surface_stations=len(stations),opposing_surface_separation_range_m=[float(gaplen.min()),float(gaplen.max())],surface_clear_points=int(safe.sum()),total_diagnostic_points=len(proposed),clearance_failure_witnesses=failures,closed_source_inside_witnesses=inside_witnesses,SP_nearest_projection_discontinuity_witnesses=jumps,
      proposed_medial_polyline_valid=False,continuous_edge_proof_not_run_because_point_gate_failed=bool(failures),source_face_adjacency_not_arterial_clearance_proof=True,scope='One literal shortest SG face-centroid chain between these exact sampled endpoints; another anatomically supported face strip remains possible. No global anatomical infeasibility claimed.',
      next_action='Review exact failed SG/SP face neighborhoods and their opposing separation. Establish whether the chosen anatomical SG/SP portal band is on the same continuous external interval as the native parent entry; no repeated grid or radius relaxation.',no_C2_or_emitted_geometry_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);print(json.dumps(dict(complete=complete,face_chain=chain,source_surface_stations=len(stations),gap_mm=(np.array(result['opposing_surface_separation_range_m'])*1000).tolist(),clearance_failure_points=len(failures),inside_sources=sorted({r['part_id']for r in inside_witnesses}),SP_projection_jumps=len(jumps),minimum_distance_mm=float(dist.min()*1000),elapsed_seconds=result['elapsed_seconds'])))

if __name__=='__main__':main()
