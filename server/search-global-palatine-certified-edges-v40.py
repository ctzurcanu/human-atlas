"""Actual-calibre graph search with lazy continuous surface-distance certificates."""
from pathlib import Path
import json,sys,time,heapq,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import AtlasReader,Mesh

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';I=D/'global-palatine-graph-v33';O=D/'global-palatine-actual-radius-v39';E=O/'field-extension-evidence.json';e=json.loads(E.read_bytes());assert e['completed']and e['total_sources']==179;C=O/'construction-contract.json';c=json.loads(C.read_bytes());P=O/'edge-search-contract-v40.json';assert not P.exists();source=ROOT/c['source_catalogue'];assert digest_file(source)==c['source_catalogue_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';B=D.parent;plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());runtime=json.loads((I/'additional-participation-frame-v37.json').read_bytes());assert digest_file(I/'additional-participation-frame-v37.json')==c['runtime_classification_sha256']
    z=np.load(I/'clearance-domain.npz');q=z['points'];shape=tuple(z['shape']);step=float(z['step']);field=ROOT/e['current_field']['path'];assert digest_file(field)==e['current_field']['sha256'];f=np.load(field);clear=f['clearance_lower_bound_m'];threshold=c['radius_m']+c['margin_m'];valid=clear>=threshold;start=np.array(c['unchanged_start_m']);end=np.array(c['unchanged_recipient']['point_m']);sgband=np.array(c['unchanged_SG_SP_band_m']);upper=np.array(c['unchanged_upper_margin_witnesses_m']);N=len(q)
    contract=dict(schema='native-global-palatine-certified-edge-search-contract/v40',source_catalogue_sha256=c['source_catalogue_sha256'],clearance_field_sha256=e['current_field']['sha256'],source_count=179,full_radius_plus_margin_m=threshold,source_station=8,unchanged_recipient=end.tolist(),ordered_sets=['station8','unchanged SG/SP band','unchanged superior-constrictor upper witnesses','unchanged PP endpoint'],
      method='A* on actual-radius nodes. Whole-edge certificates use endpoint 1-Lipschitz bounds where sufficient; uncertain edges use adaptive exact Float64 source-triangle nearest-distance lower bounds on every subinterval. No grid-half-edge amount is added to vessel radius. Parent/PP prospective interface spheres are analytically clipped and their outside intervals receive separate literal-source certificates.',
      node_occupancy_basis='Saved179-source signed/unsigned clearance field; valid closed-source occupancy establishes outside endpoints except within explicit parent/PP interface masks. A surface-clear segment cannot enter a valid closed volume without crossing its boundary.',
      distances='Persistent libigl Float64 AABB over every literal nearby triangle. Triangles outside the whole expanded domain have a conservative AABB separation greater than the required clearance.',
      known_recipient_connectors=str((I/'analytic-segment-audit-v35.json').relative_to(ROOT)),known_recipient_connectors_sha256=digest_file(I/'analytic-segment-audit-v35.json'),
      strict_parent_mask_radius_m=.003,strict_PP_recipient_mask_radius_m=.002,parent_source_fixed=True,other3_routes_literal=c['unchanged_other_route_hashes'],open_source_volume_and_complete_oral_fluid_unknown=True,
      internal_deadline_seconds=49,max_adaptive_depth=12,source_bound_progress_checkpoint=True,no_C2_or_emitted_geometry_pack_install=True)
    atomic_json(P,contract)
    spec=importlib.util.spec_from_file_location('palatine_v40_classifier',ROOT/'server/check-facial-cervical-feasibility.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h)
    metadata=e['reused_prior_sources']+e['new_sources'];meshes={};classifiers={};bindings=[];source_endpoint_inside=[]
    for r in metadata:
        pid=r['part_id']
        if pid in parts:m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=c['source_catalogue_sha256']);allow=parts[pid].get('system')in{'arterial','venous'}
        elif pid=='CANDIDATE:glandular-0:passed-fixed-envelope':
            s=plan['extra_fixed_meshes'][0];assert digest_file(ROOT/s['mesh'])==s['mesh_sha256'];m=Mesh.load(ROOT/s['mesh']);allow=True
        else:
            name=pid.split(':')[1];a=np.load(D/(name+'-fields.npz'));m=Mesh(a['fields'][0],a['faces'],frame);allow=True
        assert m.geometry_sha==r['geometry_sha256'];meshes[pid]=m;cl,mode=h.occupancy_classifier(m,allow);classifiers[pid]=cl;bindings.append(dict(part_id=pid,geometry_sha256=m.geometry_sha,classifier=mode))
        if cl is not None:
            winds=abs(igl.winding_number(cl.vertices,cl.faces,np.array([start,end])))>.99
            for k,name in enumerate(['start','recipient']):
                if winds[k]and not((name=='start'and pid=='ZA:Facial artery.l')or(name=='recipient'and pid=='ZA:Palatopharyngeus muscle.l')):source_endpoint_inside.append(dict(endpoint=name,part_id=pid,inside_verified_classifier=True))
    if source_endpoint_inside:
        atomic_json(O/'edge-search-evidence-v40.json',dict(completed=True,route_found=False,endpoint_occupied_witnesses=source_endpoint_inside,no_route_or_geometry_constructed=True,elapsed_seconds=time.monotonic()-began));print(json.dumps(dict(endpoint_occupied_witnesses=source_endpoint_inside)));return
    # Literal triangle subsets are pruned only by a complete domain AABB expanded
    # by the full radius+margin. This cannot omit a surface that blocks a query.
    domain_low=q.min(0)-threshold;domain_high=q.max(0)+threshold;vs=[];fs=[];owners=[];native_faces=[];offset=0
    parent_id='ZA:Facial artery.l';pp_id='ZA:Palatopharyngeus muscle.l'
    for source_index,r in enumerate(metadata):
        pid=r['part_id'];m=meshes[pid]
        if pid in [parent_id,pp_id]:continue
        tr=m.triangles;selected=np.flatnonzero(np.all(tr.max(1)>=domain_low,axis=1)&np.all(tr.min(1)<=domain_high,axis=1))
        if not len(selected):continue
        face=m.faces[selected];used,inv=np.unique(face,return_inverse=True);vs.append(m.vertices[used]);fs.append(inv.reshape(-1,3)+offset);offset+=len(used);owners.append(np.full(len(face),source_index,np.int32));native_faces.append(selected)
    vertices=np.concatenate(vs);faces=np.concatenate(fs);owner=np.concatenate(owners);original_face=np.concatenate(native_faces);tree=igl.AABB();tree.init(vertices,faces);trees={'composite':(tree,vertices,faces)}
    for pid in [parent_id,pp_id]:m=meshes[pid];t=igl.AABB();t.init(m.vertices,m.faces);trees[pid]=(t,m.vertices,m.faces)
    atomic_npz(O/'literal-surface-query-bank-v40.npz',vertices=vertices,faces=faces,source_index=owner,native_face=original_face,frame=np.array(frame))
    witness_queries=[];stats=dict(fast_endpoint_edges=0,adaptive_edges=0,adaptive_point_queries=0,blocked_edges=0,unresolved_edges=0);edge_cache={};edge_rows=[]
    def dist(key,point):
        t,v,ff=trees[key];d2,fi,near=t.squared_distance(v,ff,point[None]);d=float(np.sqrt(d2[0]));si=int(owner[fi[0]])if key=='composite'else next(i for i,r in enumerate(metadata)if r['part_id']==key);face_id=int(original_face[fi[0]])if key=='composite'else int(fi[0]);return d,si,face_id,near[0]
    def outside(a,b,center,r):
        u=b-a;A=float(u@u);v=a-center;B_=float(2*v@u);C_=float(v@v-r*r);disc=B_*B_-4*A*C_
        if disc<=0:return [(0.,1.)]if np.linalg.norm((a+b)/2-center)>=r else []
        roots=sorted([(-B_-np.sqrt(disc))/(2*A),(-B_+np.sqrt(disc))/(2*A)]);lo=max(0.,roots[0]);hi=min(1.,roots[1])
        if lo>=hi:return [(0.,1.)]if np.linalg.norm((a+b)/2-center)>=r else []
        return ([(0.,lo)]if lo>0 else [])+([(hi,1.)]if hi<1 else [])
    def certify(a,b,endpoint_lower=None):
        length=float(np.linalg.norm(b-a));proof=[]
        if length==0:return True,proof,'zero-length'
        if endpoint_lower is not None and endpoint_lower-length/2>=threshold:
            stats['fast_endpoint_edges']+=1;composite_intervals=[];proof.append(dict(surface='composite',method='endpoint-Lipschitz',interval=[0.,1.],lower_m=endpoint_lower-length/2))
            # Outside both prospective interface masks the saved field is a
            # lower bound on every source, including parent and PP. Thus this
            # is already a complete full-source whole-edge certificate.
            if min(np.linalg.norm(a-start),np.linalg.norm(b-start))>=.003 and min(np.linalg.norm(a-end),np.linalg.norm(b-end))>=.002:
                proof[0]['surface']='all179';return True,proof,'continuous-certified'
        else:stats['adaptive_edges']+=1;composite_intervals=[(0.,1.)]
        tasks=[('composite',a,b,composite_intervals),(parent_id,a,b,outside(a,b,start,.003)),(pp_id,a,b,outside(a,b,end,.002))]
        for key,a,b,intervals in tasks:
            for lo,hi in intervals:
                stack=[(lo,hi,0)]
                while stack:
                    lo,hi,depth=stack.pop();mid=(lo+hi)/2;point=a+mid*(b-a);d,si,face_id,surf=dist(key,point);stats['adaptive_point_queries']+=1;roundoff=512*np.finfo(float).eps*max(1.,float(np.max(abs(point))));lower=d-(hi-lo)*length/2-roundoff
                    item=dict(surface=key,method='midpoint-Lipschitz',interval=[lo,hi],sample_parameter=mid,distance_m=d,lower_m=lower,source_index=si,face=face_id,source_point_m=surf.tolist())
                    if lower>=threshold:proof.append(item);continue
                    if d<threshold:proof.append(item);stats['blocked_edges']+=1;return False,proof,'actual-distance-below-required-clearance'
                    if depth>=12:proof.append(item);stats['unresolved_edges']+=1;return False,proof,'subinterval-not-certified-at-depth-limit'
                    stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
        return True,proof,'continuous-certified'
    # Exact literal v35 endpoints are retained; only the new80 source additions
    # can invalidate a previously passing recipient segment.
    old=json.loads((I/'analytic-segment-audit-v35.json').read_bytes());goal_rows=[];goals={};alltree=cKDTree(q)
    for r in old['segments']:
        if not r['name'].startswith('recipient-')or not r['surface_and_valid_closed_source_gate_passed']:continue
        a=np.array(r['first_m']);b=np.array(r['second_m']);d,idx=alltree.query(b);assert d<1e-14 and np.array_equal(a,end)
        ok,proof,reason=certify(a,b);goal_rows.append(dict(name=r['name'],node=int(idx),passed=ok,node_passes_new179_source_field=bool(valid[idx]),reason=reason,proof=proof))
        if ok and valid[idx]:goals[int(idx)]=float(np.linalg.norm(b-a))
    starts=[];valid_ids=np.flatnonzero(valid);ds,ix=cKDTree(q[valid_ids]).query(start,k=60)
    for d,k in zip(ds,ix):
        if d>.003:continue
        idx=int(valid_ids[k]);ok,proof,reason=certify(start,q[idx]);witness_queries.append(dict(role='start',node=idx,passed=ok,reason=reason,proof=proof))
        if ok:starts.append((idx,float(d)))
    portal_set=set(np.flatnonzero(valid&(cKDTree(sgband).query(q)[0]<=.002)));upper_set=set(np.flatnonzero(valid&(cKDTree(upper).query(q)[0]<=.003)&(q[:,1]>=1.677)));portal_points=q[list(portal_set)];upper_points=q[list(upper_set)]
    if not starts or not goals or not portal_set or not upper_set:
        atomic_json(O/'edge-search-evidence-v40.json',dict(completed=True,route_found=False,source_bindings=bindings,start_connectors=witness_queries,recipient_connectors=goal_rows,set_counts=[len(starts),len(portal_set),len(upper_set),len(goals)],reason='Missing certified endpoint connector or actual-clearance milestone set',elapsed_seconds=time.monotonic()-began,no_C2_or_tube_or_install=True));print(json.dumps(dict(set_counts=[len(starts),len(portal_set),len(upper_set),len(goals)],elapsed_seconds=time.monotonic()-began)));return
    dp=cKDTree(portal_points).query(q)[0];du=cKDTree(upper_points).query(q)[0];de=np.linalg.norm(q-end,axis=1);pu=float(cKDTree(upper_points).query(portal_points)[0].min());ue=float(np.linalg.norm(upper_points-end,axis=1).min());heuristics=np.array([dp+pu+ue,du+ue,de]);distances=np.full((3,N),np.inf);previous=np.full((3,N),-1,dtype=np.int64);queue=[]
    for idx,d in starts:distances[0,idx]=d;heapq.heappush(queue,(d+heuristics[0,idx],d,idx,0))
    strides=[shape[1]*shape[2],shape[2],1];pops=0;reached=None;completed=True
    while queue:
        priority,cost,idx,stage=heapq.heappop(queue)
        if cost!=distances[stage,idx]:continue
        pops+=1
        if pops%100==0 and time.monotonic()-began>49:completed=False;break
        if stage==2 and idx in goals:reached=(stage,idx);break
        xyz=np.unravel_index(idx,shape)
        for axis in range(3):
            for direction in [-1,1]:
                if not 0<=xyz[axis]+direction<shape[axis]:continue
                other=idx+direction*strides[axis]
                if not valid[other]:continue
                pair=(min(idx,other),max(idx,other))
                if pair not in edge_cache:
                    ok,proof,reason=certify(q[pair[0]],q[pair[1]],float(min(clear[list(pair)])));edge_cache[pair]=ok
                    if proof and(reason!='continuous-certified'or len(proof)>3 or proof[0]['method']!='endpoint-Lipschitz'):edge_rows.append(dict(nodes=list(pair),passed=ok,reason=reason,proof=proof))
                if not edge_cache[pair]:continue
                ns=stage
                if stage==0 and other in portal_set:ns=1
                elif stage==1 and other in upper_set:ns=2
                nc=cost+step
                if nc<distances[ns,other]:distances[ns,other]=nc;previous[ns,other]=stage*N+idx;heapq.heappush(queue,(nc+heuristics[ns,other],nc,other,ns))
    pairs=np.array(list(edge_cache),dtype=np.int64).reshape(-1,2);statuses=np.array(list(edge_cache.values()),bool);atomic_npz(O/'edge-search-progress-v40.npz',distances=distances,previous=previous,queue=np.array(queue).reshape(-1,4),edge_pairs=pairs,edge_pass=statuses)
    states=[]
    if reached:
        st,idx=reached
        while True:
            states.append((st,idx));prior=int(previous[st,idx])
            if prior<0:break
            st,idx=divmod(prior,N)
        states=states[::-1];points=np.vstack([start,q[[i for s,i in states]],end]);atomic_npz(O/'ordered-global-route-v40.npz',points=points,states=np.array(states),frame=np.array(frame))
    result=dict(schema='native-global-palatine-certified-edge-search/v40',completed=completed,route_found=bool(reached),source_catalogue_sha256=c['source_catalogue_sha256'],contract_sha256=digest_file(P),clearance_field_sha256=e['current_field']['sha256'],query_bank_sha256=digest_file(O/'literal-surface-query-bank-v40.npz'),source_bindings=bindings,set_counts=[len(starts),len(portal_set),len(upper_set),len(goals)],start_connectors=witness_queries,recipient_connectors=goal_rows,visited_states=pops,highest_reached_milestone=max(np.flatnonzero(np.isfinite(distances).any(1)),default=-1),edge_count=len(edge_cache),statistics=stats,nontrivial_edge_certificates=edge_rows,progress_sha256=digest_file(O/'edge-search-progress-v40.npz'),route_points=len(states)+2 if reached else 0,
      route_sha256=digest_file(O/'ordered-global-route-v40.npz')if reached else None,remaining_queue=len(queue),reason='route found or certified reachable frontier exhausted'if completed else 'bounded search checkpoint; not an anatomical infeasibility result',elapsed_seconds=time.monotonic()-began,no_C2_or_emitted_geometry_pack_install=True,open_physical_volumes_and_complete_oral_fluid_unknown=True)
    atomic_json(O/'edge-search-evidence-v40.json',result);assert digest_file(source)==c['source_catalogue_sha256'];print(json.dumps({k:result[k]for k in ['completed','route_found','set_counts','visited_states','highest_reached_milestone','edge_count','statistics','remaining_queue','elapsed_seconds']}))

if __name__=='__main__':main()
