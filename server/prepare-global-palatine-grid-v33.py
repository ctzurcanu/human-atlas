"""One saved, native global graph proposal. This is not an emitted tube or anatomy acceptance."""
from pathlib import Path
import sys, json, time, importlib.util, heapq
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader

def main():
    began=time.monotonic();limit=49.
    B=ROOT/'.local-models/geometry-workbench/circulatory-completion';D=B/'male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'global-palatine-graph-v33';assert not O.exists();O.mkdir()
    expected='d73c260cf6ed0752af5b3f420e97b11f4f180bea100aaf2f70f9082bf9a675de'
    binding=json.loads((D/'source-binding-v19/evidence.json').read_bytes());source=ROOT/binding['snapshot'];live=ROOT/binding['current_catalogue'];assert digest_file(source)==expected==digest_file(live)
    atlas=json.loads(source.read_bytes());parts={r['id']:r for r in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres'
    plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());proof=json.loads((B/'male-left-facial-cervical-routes-palatine-active-c2-v28/literal-source-rebase.json').read_bytes())
    parent=Mesh.load(ROOT/plan['fixed_parent']['mesh']);rings=np.array(plan['fixed_parent']['source_station_components'][0]['rings']);start=parent.vertices[rings[8]].mean(0);route=plan['routes'][0];end=np.array(route['recipient']['point_m']);radius=.0005501;margin=.00007;step=.001
    portals=json.loads((D/'global-corridor-analysis-v31/evidence.json').read_bytes())['coarse_portal_landmarks']
    # A band of actual opposing SG/SP witnesses, not a list of pinned guide vertices.
    band=np.array([r['point_m']for r in portals if 1.640<=r['point_m'][1]<=1.654])
    spec=importlib.util.spec_from_file_location('global_palatine_classifier_v33',ROOT/'server/check-facial-cervical-feasibility.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    fixed={};source_rows=[]
    for row in proof['fixed_blocker_geometry']:
        pid=row['part_id']
        if pid in parts:m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=expected);allow=parts[pid].get('system')in{'arterial','venous'}
        else:
            extra=next(r for r in plan['extra_fixed_meshes']if r['part_id']==pid);assert digest_file(ROOT/extra['mesh'])==extra['mesh_sha256'];m=Mesh.load(ROOT/extra['mesh']);allow=True
        assert m.geometry_sha==row['geometry_sha256'];c,mode=h.occupancy_classifier(m,allow);fixed[pid]=(m,c,mode)
    for pid in ['ZA:Gingiva','ZA:Soft palate','O3M:Trapezius muscle.l']:
        m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=expected);c,mode=h.occupancy_classifier(m,False);fixed[pid]=(m,c,mode)
    assert len(fixed)==96
    for name in ['tonsillar','glandular-1']:
        z=np.load(D/(name+'-fields.npz'));m=Mesh(z['fields'][0],z['faces'],frame);c,mode=h.occupancy_classifier(m,True);assert c is not None;fixed['CANDIDATE:'+name+':fixed-envelope']=(m,c,mode)
    c,mode=h.occupancy_classifier(parent,True);assert c is not None;fixed['ZA:Facial artery.l']=(parent,c,mode)
    sc=fixed['ZA:Superior pharyngeal constrictor.l'][0]
    # The upper witness set is derived from literal source vertices within each
    # 2 mm lateral bin. It is a coarse margin location, not an anatomical edge.
    upper=[]
    for x in np.arange(.004,.0301,.002):
        ids=np.flatnonzero(abs(sc.vertices[:,0]-x)<=.001)
        if len(ids):upper.append(sc.vertices[ids[np.argmax(sc.vertices[ids,1])]])
    upper=np.array(upper)
    low=np.minimum.reduce([start,end,band.min(0),upper.min(0)])-np.array([.006,.004,.010]);high=np.maximum.reduce([start,end,band.max(0),upper.max(0)])+np.array([.009,.004,.009])
    axes=[np.arange(low[k],high[k]+step/2,step)for k in range(3)];shape=tuple(len(a)for a in axes);q=np.stack(np.meshgrid(*axes,indexing='ij'),-1).reshape(-1,3)
    contract=dict(schema='native-global-palatine-graph-contract/v33',source_catalogue=str(source.relative_to(ROOT)),source_catalogue_sha256=expected,registry_revision='31fab36e6897adaec9c8ddca201261607902203d75c64f3b214e57d329f460c2',
       authority_inventory= 'reports/ta98-male/facial-arterial-family-contract-v1.json',authority_inventory_sha256=digest_file(ROOT/'reports/ta98-male/facial-arterial-family-contract-v1.json'),
       term='A12.2.05.021',parent_term='A12.2.05.020',representation='future finite vascular wall plus separate lumen, not opening plugs',stage='global free-space graph proposal only',
       fixed_parent_id='ZA:Facial artery.l',fixed_native_station=8,start_m=start.tolist(),recipient=route['recipient'],source_profile_radius_m=radius,minimum_clearance_m=margin,
       parent_interface_mask='Only centers inside 3 mm of unchanged station8 may overlap the native parent; this is a declared prospective shared junction, not a completed or accepted interface.',
       recipient_mask='Only centers inside 2 mm of unchanged PP endpoint may occupy PP tissue. No oral fluid space is constructed or authorized.',
       domain_bounds_m=[low.tolist(),high.tolist()],spacing_m=step,shape=shape,grid_points=len(q),ordered_milestones=['source station8','SG/SP opposing-surface band','literal superior-constrictor upper-margin witness band','unchanged PP recipient'],
       SG_SP_band_m=band.tolist(),superior_constrictor_upper_witnesses_m=upper.tolist(),unmodified_other_routes={name:digest_file(D/(name+'-fields.npz'))for name in ['tonsillar','glandular-0','glandular-1']},
       algorithm='6-neighbor uniform grid Dijkstra through ordered anatomical sets; intermediary nodes are free. Every axis edge uses min endpoint clearance - h/2 as a 1-Lipschitz surface-distance lower bound, so an accepted edge retains the full source radius plus margin. This is a graph, not a C2 tube.',
       complete_96_declared_default_sources=True,open_physical_occupancy_unknowns_remain=True,complete_oral_fluid_exclusion_unresolved=True,not_accepted_anatomy=True,no_pack_or_install=True)
    atomic_json(O/'construction-contract.json',contract)
    # Rebind actual instanced hair to this domain, never the carrier hull.
    hp=parts['ZA:Scalp hair:metric'];field=hp['hairField'];hp_path=reader.asset(field['url']);assert digest_file(hp_path)==field['sha256'];raw=np.memmap(hp_path,dtype='u1',mode='r');roots=reader.view(raw,field['positions'],field['count'],'<f4',3).astype(float);packed=reader.view(raw,field['normals'],field['count'],'<i2',3);dirs=np.clip(packed.astype(float)/32767.,-1,1);dirs/=np.linalg.norm(dirs,axis=1)[:,None];tips=roots+dirs*field['length'];r=field['diameter']/2;eps=np.finfo(np.float32).eps;halo=float(4*(64*eps/(1-64*eps))*(np.max(abs(roots))+field['length']+r));extent=r*np.sqrt(np.maximum(0,1-dirs**2));hl=np.minimum(roots,tips)-extent-halo;hh=np.maximum(roots,tips)+extent+halo;hair_ids=np.flatnonzero(np.all(hh>=low-radius,axis=1)&np.all(hl<=high+radius,axis=1));atomic_npz(O/'hair-domain-overlap.npz',shaft_ids=hair_ids)
    hair=dict(sha256=field['sha256'],count=int(field['count']),diameter_m=field['diameter'],length_m=field['length'],float32_halo_m=halo,domain_shaft_aabb_overlap_count=int(len(hair_ids)),no_carrier_used=True)
    # A field of lower bounds is enough for graph certification. Distant-source
    # AABBs are lower bounds, queried regions get exact source distances.
    clearance=np.full(len(q),.010);limiting=np.full(len(q),-1,dtype=np.int32);complete=True
    for n,(pid,(m,c,mode))in enumerate(fixed.items()):
        if time.monotonic()-began>limit:complete=False;break
        box=np.linalg.norm(np.maximum(m.bounds[0]-q,0)+np.maximum(q-m.bounds[1],0),axis=1);take=np.flatnonzero(box<.0022);values=np.minimum(box,.010)
        inside=np.zeros(len(take),bool)
        if len(take):
            d2,faces,near=igl.point_mesh_squared_distance(q[take],m.vertices,m.faces);dist=np.sqrt(d2)
            if c is not None:inside=abs(igl.winding_number(c.vertices,c.faces,q[take]))>.99;dist[inside]*=-1
            values[take]=dist
        if pid=='ZA:Palatopharyngeus muscle.l':values[np.linalg.norm(q-end,axis=1)<=.002]=.010
        if pid=='ZA:Facial artery.l':values[np.linalg.norm(q-start,axis=1)<=.003]=.010
        new=values<clearance;clearance[new]=values[new];limiting[new]=n
        source_rows.append(dict(part_id=pid,geometry_sha256=m.geometry_sha,classifier=mode,queried=len(take),closed_source_inside_grid_count=int(inside.sum()),radius_and_margin_violations=int((values<radius+margin).sum())))
    # Unclassified newly retrieved sources are retained in the evidence, not
    # quietly considered empty. A proposal cannot establish their clearance.
    additional=[]
    for p in atlas['parts']:
        if p['id']in fixed or p['id']=='ZA:Scalp hair:metric' or not p.get('bounds') or p.get('suppressed'):continue
        bb=np.array(p['bounds'])
        if np.all(bb[1]>=low-radius)&np.all(bb[0]<=high+radius):additional.append(dict(part_id=p['id'],system=p.get('system'),defaultHidden=p.get('defaultHidden',False),suppressed=p.get('suppressed',False),bounds_m=bb.tolist()))
    atomic_npz(O/'clearance-domain.npz',points=q,clearance_lower_bound_m=clearance,limiting_source=limiting,shape=np.array(shape),low=low,step=np.array(step),frame=np.array(frame))
    grid_output=dict(schema='native-global-palatine-graph/v33',source_catalogue_sha256=expected,contract_sha256=digest_file(O/'construction-contract.json'),domain_sha256=digest_file(O/'clearance-domain.npz'),fixed_source_count=len(fixed),source_queries=source_rows,complete_declared_context=complete,hair=hair,newly_retrieved_unclassified_sources=additional,elapsed_seconds=time.monotonic()-began,not_a_C2_guide=True,no_emitted_geometry=True,no_pack_or_install=True,continuous_graph_edge_certificate_method=contract['algorithm'],physical_volume_unknowns=[x['part_id']for x in source_rows if x['classifier']in ['source_is_not_a_closed_occupied_volume','invalid_native_conduit_topology','unresolved_native_conduit_boundary','temporary_conduit_classifier_is_not_closed']],complete_oral_fluid_exclusion_unresolved=True)
    atomic_json(O/'domain-evidence.json',grid_output)
    if not complete or len(hair_ids):
        atomic_json(O/'graph-evidence.json',dict(completed=False,reason='domain construction internal deadline'if not complete else 'actual hair shaft domain overlap requires explicit tests',elapsed_seconds=time.monotonic()-began));print(json.dumps(dict(completed=False,queried_sources=len(source_rows),elapsed_seconds=time.monotonic()-began)));return
    # Certify whole axis edges from endpoint lower bounds. Nearest start/goal
    # connector segments use the same lower-bound inequality and all 99 sources.
    valid=clearance>=radius+margin+step/2;coords=np.arange(len(q)).reshape(shape);valid_grid=valid.reshape(shape);qtree=cKDTree(q[valid]);valid_ids=np.flatnonzero(valid)
    def endpoint_candidates(point,exemption):
        ds,ix=qtree.query(point,k=min(100,len(valid_ids)));ids=valid_ids[np.atleast_1d(ix)];ds=np.atleast_1d(ds);out=[]
        values=[]
        for pid,(m,c,mode)in fixed.items():
            if pid==exemption:continue
            box=float(np.linalg.norm(np.maximum(m.bounds[0]-point,0)+np.maximum(point-m.bounds[1],0)))
            if box>=.004:values.append(box);continue
            dist=float(np.sqrt(igl.point_mesh_squared_distance(point[None],m.vertices,m.faces)[0][0]));
            if c is not None and abs(igl.winding_number(c.vertices,c.faces,point[None])[0])>.99:dist=-dist
            values.append(dist)
        endpoint_clear=min(values)
        for idx,d in zip(ids,ds):
            if d>(.003 if exemption=='ZA:Facial artery.l' else .002):continue
            # Sample exactly at start/goal and use the nearest-endpoint bound.
            if min(endpoint_clear,clearance[idx])-d/2>=radius+margin:out.append((int(idx),float(d)))
        return out
    starts=endpoint_candidates(start,'ZA:Facial artery.l');ends=endpoint_candidates(end,'ZA:Palatopharyngeus muscle.l')
    portal_set=set(np.flatnonzero(valid&(cKDTree(band).query(q)[0]<=.002)))
    upper_set=set(np.flatnonzero(valid&(cKDTree(upper).query(q)[0]<=.003)&(q[:,1]>=1.677)))
    goals={i:d for i,d in ends};goal_set=set(goals);distance={};prev={};queue=[];reached=None;pops=0
    for idx,cost in starts:distance[(idx,0)]=cost;heapq.heappush(queue,(cost,idx,0))
    strides=np.array([shape[1]*shape[2],shape[2],1]);milestone_counts=[len(portal_set),len(upper_set),len(ends)]
    graph_finished=True
    while queue:
        cost,idx,stage=heapq.heappop(queue);state=(idx,stage)
        if cost!=distance.get(state):continue
        pops+=1
        if pops%5000==0 and time.monotonic()-began>limit:graph_finished=False;break
        if stage==2 and idx in goal_set:reached=state;break
        xyz=np.array(np.unravel_index(idx,shape))
        for k in range(3):
            for direction in [-1,1]:
                if not 0<=xyz[k]+direction<shape[k]:continue
                neighbor=idx+direction*int(strides[k])
                if not valid[neighbor]:continue
                ns=stage
                if stage==0 and neighbor in portal_set:ns=1
                elif stage==1 and neighbor in upper_set:ns=2
                key=(neighbor,ns);nc=cost+step
                if nc<distance.get(key,np.inf):distance[key]=nc;prev[key]=state;heapq.heappush(queue,(nc,neighbor,ns))
    states=[]
    if reached:
        key=reached
        while key in prev:states.append(key);key=prev[key]
        states.append(key);states=states[::-1];points=np.vstack([start,q[[k[0]for k in states]],end]);atomic_npz(O/'ordered-global-graph-route.npz',points=points,states=np.array(states),frame=np.array(frame))
    else:points=np.empty((0,3))
    source_point_rows=[]
    for pid,(m,c,mode)in fixed.items():
        d2,f,p=igl.point_mesh_squared_distance(start[None],m.vertices,m.faces);source_point_rows.append(dict(part_id=pid,distance_m=float(np.sqrt(d2[0])),source_face=int(f[0]),source_point_m=p[0].tolist(),classifier=mode))
    result=dict(completed=graph_finished,reason='graph frontier exhausted or route found'if graph_finished else 'graph internal wall-clock deadline',route_found=bool(reached),start_connectors=len(starts),end_connectors=len(ends),free_grid_nodes=int(valid.sum()),milestone_set_counts=milestone_counts,visited_states=pops,reached_milestone=max((s for _,s in distance),default=-1),graph_length_m=float(np.linalg.norm(np.diff(points,axis=0),axis=1).sum())if reached else None,
        source_entry_closest_witnesses=sorted(source_point_rows,key=lambda x:x['distance_m'])[:12],ordered_route_points=int(len(points)),route_file_sha256=digest_file(O/'ordered-global-graph-route.npz')if reached else None,
        clearance_scope='Complete declared96 default sources, other3 fixed branch envelope set and native parent outside prospective interface mask; additional source retrieval remains unclassified.',
        continuous_edge_clearance_lower_bound_m=radius+margin,not_a_C2_curve_or_emitted_tube=True,no_anatomical_acceptance=True,no_pack_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'graph-evidence.json',result);assert digest_file(source)==expected==digest_file(live)
    print(json.dumps({k:result[k]for k in ['route_found','start_connectors','end_connectors','free_grid_nodes','milestone_set_counts','visited_states','reached_milestone','graph_length_m','elapsed_seconds']}))

if __name__=='__main__':main()
