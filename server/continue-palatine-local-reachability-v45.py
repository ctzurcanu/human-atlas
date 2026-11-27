"""Continue exactly v44's local graph for reachability, without shortest-path work."""
from pathlib import Path
import sys,json,time
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from collections import deque
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';I=D/'palatine-surface-portal-graph-v44';G=D/'global-palatine-actual-radius-v39';O=D/'palatine-surface-portal-reachability-v45';assert not O.exists();O.mkdir();E=I/'evidence.json';e=json.loads(E.read_bytes());assert not e['completed'];C=I/'construction-contract.json';c=json.loads(C.read_bytes());assert digest_file(ROOT/'.local-models/ta98-review/bowel-wall-local-atlas/male/atlas.json')==c['source_catalogue_sha256'];assert digest_file(I/'local-graph-progress.npz')==e['progress_sha256'];p=np.load(I/'local-graph-progress.npz');points=p['points'];clear=p['clearance_m'];n=len(points);assert n==704;threshold=c['radius_m']+c['margin_m'];seed_rows={int(r['local_node']):(j,r)for j,r in enumerate(e['certificates'])if r.get('grid_node')is not None and r['passed']};assert len(seed_rows)==684;cache={tuple(map(int,x)):bool(y)for x,y in zip(p['edges'],p['edge_pass'])};oldproof={tuple(r['local_nodes']):j for j,r in enumerate(e['certificates'])if r.get('local_nodes')is not None};band=np.array(c['unchanged_SG_SP_band']);target=set(map(int,np.flatnonzero(cKDTree(band).query(points)[0]<=.002)));tree=cKDTree(points);neighbors=tree.query_ball_point(points,.0015);reachable=np.zeros(n,bool);reachable[list(seed_rows)]=True;parent=np.full(n,-1,dtype=np.int64);queue=deque(sorted(seed_rows));proofrefs=[];newproofs=[];stats=dict(reused_passed=0,reused_blocked=0,new_certifications=0,adaptive_queries=0,failed=0,unresolved=0)
    contract=dict(schema='palatine-fixed-local-graph-reachability-contract/v45',source_catalogue_sha256=c['source_catalogue_sha256'],original_contract_sha256=digest_file(C),prior_incomplete_evidence_sha256=digest_file(E),prior_progress_sha256=e['progress_sha256'],fixed_local_nodes=n,certified_grid_seed_nodes=len(seed_rows),initially_unreached_nodes=n-len(seed_rows),required_SG_SP_band_nodes=len(target),radius_m=c['radius_m'],margin_m=c['margin_m'],method='Breadth-first reachability from every literal previously certified grid seed. Reuse cached certificates; only edges into previously unreached nodes can require new continuous source checks. No path optimization, new sampling or field regeneration.',source_station8_and_tissues_fixed=True,no_C2_emitted_geometry_pack_install=True,internal_work_deadline_seconds=48)
    atomic_json(O/'construction-contract.json',contract)
    bank=np.load(G/'literal-surface-query-bank-v40.npz');v40=json.loads((G/'edge-search-evidence-v40.json').read_bytes());assert digest_file(G/'literal-surface-query-bank-v40.npz')==v40['query_bank_sha256'];bvh=igl.AABB();bvh.init(bank['vertices'],bank['faces']);trees=[('composite',bvh,bank['vertices'],bank['faces'])]
    for pid,name in [('ZA:Facial artery.l','facial-left.npz'),('ZA:Palatopharyngeus muscle.l','palatopharyngeus-left.npz')]:
        m=Mesh.load(D.parent/'male-left-facial-cervical-routes-anatomical-v5'/name);assert m.geometry_sha==next(r['geometry_sha256']for r in v40['source_bindings']if r['part_id']==pid);t=igl.AABB();t.init(m.vertices,m.faces);trees.append((pid,t,m.vertices,m.faces))
    def certify(a,b,lower):
        stats['new_certifications']+=1;length=float(np.linalg.norm(b-a))
        if lower-length/2>=threshold:return True,dict(method='endpoint-Lipschitz',lower_m=lower-length/2)
        proof=[]
        for pid,t,v,f in trees:
            stack=[(0.,1.,0)]
            while stack:
                lo,hi,depth=stack.pop();mid=(lo+hi)/2;q=a+mid*(b-a);d2,face,near=t.squared_distance(v,f,q[None]);distance=float(np.sqrt(d2[0]));bound=distance-(hi-lo)*length/2-3e-13;stats['adaptive_queries']+=1
                if bound>=threshold:proof.append(dict(surface=pid,interval=[lo,hi],distance_m=distance,lower_m=bound,face=int(face[0])));continue
                if distance<threshold:stats['failed']+=1;return False,dict(method='actual-clearance-deficit',surface=pid,parameter=mid,distance_m=distance,face=int(face[0]),source_point_m=near[0].tolist())
                if depth==12:stats['unresolved']+=1;return False,dict(method='adaptive-depth-limit',surface=pid,parameter=mid,distance_m=distance,lower_m=bound,face=int(face[0]))
                stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
        return True,dict(method='adaptive-continuous-Lipschitz',proof=proof)
    found=next(iter(target&set(seed_rows)),None);complete=True;visited=0
    while queue and found is None:
        if time.monotonic()-began>48:complete=False;break
        a=queue.popleft();visited+=1
        for b in neighbors[a]:
            if reachable[b]:continue
            pair=(min(a,b),max(a,b))
            if pair in cache:
                ok=cache[pair];stats['reused_passed'if ok else'reused_blocked']+=1;ref=dict(prior_evidence_certificate_index=oldproof[pair],prior_evidence_sha256=contract['prior_incomplete_evidence_sha256'])
            else:
                ok,proof=certify(points[pair[0]],points[pair[1]],float(min(clear[list(pair)])));cache[pair]=ok;ref=dict(new_proof_index=len(newproofs));newproofs.append(dict(nodes=list(pair),passed=ok,proof=proof))
            if ok:
                reachable[b]=True;parent[b]=a;queue.append(b);proofrefs.append(dict(from_node=a,to_node=b,certificate=ref))
                if b in target:found=b;break
    atomic_npz(O/'reachability-progress.npz',points=points,reachable=reachable,parent=parent,queue=np.array(queue,dtype=np.int64),edge_pairs=np.array(list(cache),dtype=np.int64).reshape(-1,2),edge_pass=np.array(list(cache.values()),bool))
    connection=None;seed=None
    if found is not None:
        connection=[];k=found
        while True:
            connection.append(k)
            if k in seed_rows:seed=k;break
            k=int(parent[k]);assert k>=0
        connection=connection[::-1];seed_certificate=seed_rows[seed];grid_id=seed_certificate[1]['grid_node'];grid=np.load(D/'global-palatine-graph-v33/clearance-domain.npz')['points'];atomic_npz(O/'source-grid-to-SG_SP-band-connection.npz',grid_source_node=np.array(grid_id),local_node_ids=np.array(connection),points=np.vstack([grid[grid_id],points[connection]]),frame=np.array('model:local-ta98:assembled-metres'))
    result=dict(schema='palatine-fixed-local-reachability-continuation/v45',completed=complete,source_catalogue_sha256=c['source_catalogue_sha256'],contract_sha256=digest_file(O/'construction-contract.json'),prior_evidence_sha256=contract['prior_incomplete_evidence_sha256'],source_to_SG_SP_band_connection_found=found is not None,initial_seed_nodes=len(seed_rows),reached_nodes=int(reachable.sum()),visited_nodes=visited,remaining_queue=len(queue),unreached_nodes=np.flatnonzero(~reachable).tolist(),target_nodes=sorted(target),statistics=stats,used_proof_references=proofrefs,new_edge_proofs=newproofs,
      seed_certificate_index=seed_rows[seed][0]if seed is not None else None,local_path_node_ids=connection,connection_sha256=digest_file(O/'source-grid-to-SG_SP-band-connection.npz')if found is not None else None,progress_sha256=digest_file(O/'reachability-progress.npz'),outcome='Certified source-grid to required SG/SP band connection found; whole ordered path and C2 remain absent'if found is not None else 'Finite fixed surface-sample graph has no certified source-to-band connection'if complete else 'Explicit bounded continuation checkpoint; work incomplete',no_new_sampling_field_geometry_or_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);print(json.dumps({k:result[k]for k in ['completed','source_to_SG_SP_band_connection_found','initial_seed_nodes','reached_nodes','visited_nodes','remaining_queue','target_nodes','statistics','outcome','elapsed_seconds']}))

if __name__=='__main__':main()
