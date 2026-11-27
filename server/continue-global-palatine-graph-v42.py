"""Continue v40's saved reachable graph using only v41's four new certificates."""
from pathlib import Path
import json,sys,time,heapq
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from geometry_tools.model import Mesh
from geometry_tools.common import atomic_json,atomic_npz,digest_file

def outside(a,b,center,r):
    u=b-a;A=float(u@u);v=a-center;B=float(2*v@u);C=float(v@v-r*r);disc=B*B-4*A*C
    if disc<=0:return [(0.,1.)]if np.linalg.norm((a+b)/2-center)>=r else []
    roots=sorted([(-B-np.sqrt(disc))/(2*A),(-B+np.sqrt(disc))/(2*A)]);lo=max(0.,roots[0]);hi=min(1.,roots[1])
    if lo>=hi:return [(0.,1.)]if np.linalg.norm((a+b)/2-center)>=r else []
    return ([(0.,lo)]if lo>0 else [])+([(hi,1.)]if hi<1 else [])

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';O=D/'global-palatine-actual-radius-v39';E=O/'edge-search-evidence-v40.json';e=json.loads(E.read_bytes());R=O/'uncertain-edge-resolutions-v41.json';r=json.loads(R.read_bytes());c=json.loads((O/'construction-contract.json').read_bytes());assert e['completed']and not e['route_found']and not e['remaining_queue'];assert digest_file(E)==r['prior_search_sha256'];assert digest_file(ROOT/c['source_catalogue'])==c['source_catalogue_sha256'];assert digest_file(O/'edge-search-progress-v40.npz')==e['progress_sha256'];assert digest_file(O/'literal-surface-query-bank-v40.npz')==e['query_bank_sha256'];P=O/'continuation-contract-v42.json';assert not P.exists()
    atomic_json(P,dict(schema='global-palatine-saved-graph-continuation-contract/v42',source_catalogue_sha256=c['source_catalogue_sha256'],prior_terminal_search_sha256=digest_file(E),prior_progress_sha256=e['progress_sha256'],new_edge_certificates_sha256=digest_file(R),literal_query_bank_sha256=e['query_bank_sha256'],new_valid_edges=[x['nodes']for x in r['rows']if x['passed']],retained_failed_edges=[x['nodes']for x in r['rows']if not x['passed']],source_count=179,radius_m=c['radius_m'],margin_m=c['margin_m'],no_field_regeneration=True,unchanged_SG_SP_band=c['unchanged_SG_SP_band_m'],unchanged_upper_witnesses=c['unchanged_upper_margin_witnesses_m'],source_station8_fixed=True,recipient_unchanged=True,other_routes_literal=c['unchanged_other_route_hashes'],method='Resume the exact saved finite shortest-path labels and edges. Add only four analytically validated edges, seed their newly relaxed endpoints; continue certified traversal. No clearance relaxation or graph reconstruction.',no_C2_or_emitted_geometry_or_install=True))
    z=np.load(D/'global-palatine-graph-v33/clearance-domain.npz');q=z['points'];shape=tuple(z['shape']);step=float(z['step']);N=len(q);field=json.loads((O/'field-extension-evidence.json').read_bytes())['current_field'];assert digest_file(ROOT/field['path'])==field['sha256'];f=np.load(ROOT/field['path']);clear=f['clearance_lower_bound_m'];threshold=c['radius_m']+c['margin_m'];valid=clear>=threshold;start=np.array(c['unchanged_start_m']);end=np.array(c['unchanged_recipient']['point_m']);bank=np.load(O/'literal-surface-query-bank-v40.npz');tree=igl.AABB();tree.init(bank['vertices'],bank['faces']);trees={'composite':(tree,bank['vertices'],bank['faces'])}
    parent_id='ZA:Facial artery.l';pp_id='ZA:Palatopharyngeus muscle.l'
    for pid,name in [(parent_id,'facial-left.npz'),(pp_id,'palatopharyngeus-left.npz')]:
        m=Mesh.load(D.parent/'male-left-facial-cervical-routes-anatomical-v5'/name);binding=next(x for x in e['source_bindings']if x['part_id']==pid);assert m.geometry_sha==binding['geometry_sha256'];t=igl.AABB();t.init(m.vertices,m.faces);trees[pid]=(t,m.vertices,m.faces)
    progress=np.load(O/'edge-search-progress-v40.npz');distance=progress['distances'].copy();prev=progress['previous'].copy();cache={tuple(map(int,p)):bool(v)for p,v in zip(progress['edge_pairs'],progress['edge_pass'])};queue=[];edge_rows=[];stats=dict(fast_endpoint_edges=0,adaptive_edges=0,adaptive_point_queries=0,blocked_edges=0,unresolved_edges=0)
    portal_set=set(np.flatnonzero(valid&(cKDTree(np.array(c['unchanged_SG_SP_band_m'])).query(q)[0]<=.002)));upper_set=set(np.flatnonzero(valid&(cKDTree(np.array(c['unchanged_upper_margin_witnesses_m'])).query(q)[0]<=.003)&(q[:,1]>=1.677)));goals={x['node']:float(np.linalg.norm(q[x['node']]-end))for x in e['recipient_connectors']if x['passed']and x['node_passes_new179_source_field']}
    def dist(key,point):
        t,v,ff=trees[key];d2,fi,near=t.squared_distance(v,ff,point[None]);index=int(fi[0]);si=int(bank['source_index'][index])if key=='composite'else next(i for i,x in enumerate(e['source_bindings'])if x['part_id']==key);face=int(bank['native_face'][index])if key=='composite'else index;return float(np.sqrt(d2[0])),si,face,near[0]
    def certify(a,b,endpoint_lower):
        length=float(np.linalg.norm(b-a));proof=[]
        if endpoint_lower-length/2>=threshold:
            stats['fast_endpoint_edges']+=1;intervals=[];proof.append(dict(surface='composite',method='endpoint-Lipschitz',interval=[0.,1.],lower_m=endpoint_lower-length/2))
            if min(np.linalg.norm(a-start),np.linalg.norm(b-start))>=.003 and min(np.linalg.norm(a-end),np.linalg.norm(b-end))>=.002:return True,proof,'continuous-certified'
        else:stats['adaptive_edges']+=1;intervals=[(0.,1.)]
        tasks=[('composite',intervals),(parent_id,outside(a,b,start,.003)),(pp_id,outside(a,b,end,.002))]
        for key,intervals in tasks:
            for lo,hi in intervals:
                stack=[(lo,hi,0)]
                while stack:
                    lo,hi,depth=stack.pop();mid=(lo+hi)/2;point=a+mid*(b-a);d,si,face,surf=dist(key,point);stats['adaptive_point_queries']+=1;roundoff=512*np.finfo(float).eps*max(1.,float(np.max(abs(point))));lower=d-(hi-lo)*length/2-roundoff;item=dict(surface=key,method='midpoint-Lipschitz',interval=[lo,hi],sample_parameter=mid,distance_m=d,lower_m=lower,source_index=si,face=face,source_point_m=surf.tolist())
                    if lower>=threshold:proof.append(item);continue
                    if d<threshold:proof.append(item);stats['blocked_edges']+=1;return False,proof,'actual-distance-below-required-clearance'
                    if depth>=12:proof.append(item);stats['unresolved_edges']+=1;return False,proof,'subinterval-not-certified-at-depth-limit'
                    stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
        return True,proof,'continuous-certified'
    def nextstage(stage,idx):
        if stage==0 and idx in portal_set:return 1
        if stage==1 and idx in upper_set:return 2
        return stage
    seeded=[]
    for row in r['rows']:
        a,b=map(int,row['nodes']);pair=(a,b);assert pair in cache and not cache[pair];cache[pair]=row['passed']
        if not row['passed']:continue
        for x,y in [(a,b),(b,a)]:
            for stage in range(3):
                old=distance[stage,x]
                if not np.isfinite(old):continue
                ns=nextstage(stage,y);cost=old+step
                if cost<distance[ns,y]:distance[ns,y]=cost;prev[ns,y]=stage*N+x;heapq.heappush(queue,(cost,y,ns));seeded.append(dict(from_node=x,to_node=y,stage=ns,cost=cost))
    strides=[shape[1]*shape[2],shape[2],1];pops=0;reached=None;complete=True
    while queue:
        cost,idx,stage=heapq.heappop(queue)
        if cost!=distance[stage,idx]:continue
        pops+=1
        if pops%100==0 and time.monotonic()-began>49:complete=False;break
        if stage==2 and idx in goals:reached=(stage,idx);break
        xyz=np.unravel_index(idx,shape)
        for axis in range(3):
            for direction in [-1,1]:
                if not 0<=xyz[axis]+direction<shape[axis]:continue
                other=idx+direction*strides[axis]
                if not valid[other]:continue
                pair=(min(idx,other),max(idx,other))
                if pair not in cache:
                    ok,proof,reason=certify(q[pair[0]],q[pair[1]],float(min(clear[list(pair)])));cache[pair]=ok;edge_rows.append(dict(nodes=list(pair),passed=ok,reason=reason,proof=proof))
                if not cache[pair]:continue
                ns=nextstage(stage,other);nc=cost+step
                if nc<distance[ns,other]:distance[ns,other]=nc;prev[ns,other]=stage*N+idx;heapq.heappush(queue,(nc,other,ns))
    atomic_npz(O/'continued-progress-v42.npz',distances=distance,previous=prev,queue=np.array(queue).reshape(-1,3),edge_pairs=np.array(list(cache),dtype=np.int64),edge_pass=np.array(list(cache.values()),bool))
    if reached:
        states=[];st,idx=reached
        while True:
            states.append((st,idx));prior=int(prev[st,idx])
            if prior<0:break
            st,idx=divmod(prior,N)
        states=states[::-1];points=np.vstack([start,q[[i for s,i in states]],end]);atomic_npz(O/'ordered-global-route-v42.npz',points=points,states=np.array(states),frame=np.array('model:local-ta98:assembled-metres'))
    result=dict(schema='global-palatine-saved-certified-graph-continuation/v42',completed=complete,route_found=bool(reached),source_catalogue_sha256=c['source_catalogue_sha256'],contract_sha256=digest_file(P),prior_search_sha256=digest_file(E),new_analytic_certificates_sha256=digest_file(R),prior_edges=len(progress['edge_pass']),newly_seeded_labels=seeded,new_visited_states=pops,total_edges=len(cache),new_edge_certificates=edge_rows,statistics=stats,remaining_queue=len(queue),highest_reached_milestone=int(max(np.flatnonzero(np.isfinite(distance).any(1)),default=-1)),progress_sha256=digest_file(O/'continued-progress-v42.npz'),route_sha256=digest_file(O/'ordered-global-route-v42.npz')if reached else None,
      outcome='Valid certified discrete route found; no C2/emitted tube acceptance'if reached else 'Certified reachable discrete frontier exhausted; not a proof of global anatomical infeasibility'if complete else 'Saved bounded continuation checkpoint',no_field_regeneration=True,no_C2_or_emitted_geometry_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'continued-search-evidence-v42.json',result);print(json.dumps({k:result[k]for k in ['completed','route_found','new_visited_states','total_edges','statistics','remaining_queue','highest_reached_milestone','outcome','elapsed_seconds']}))

if __name__=='__main__':main()
