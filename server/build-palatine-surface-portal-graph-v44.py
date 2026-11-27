"""One source-surface adaptive SG/SP portal connection stage, with no new full field."""
from pathlib import Path
import json,sys,time,heapq,importlib.util
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import numpy as np
import igl
from scipy.spatial import cKDTree
from geometry_tools.common import atomic_json,atomic_npz,digest_file
from geometry_tools.model import Mesh,AtlasReader

def main():
    began=time.monotonic();D=ROOT/'.local-models/geometry-workbench/circulatory-completion/male-left-facial-cervical-finite-network-cervical-c2-v17';I=D/'global-palatine-graph-v33';G=D/'global-palatine-actual-radius-v39';O=D/'palatine-surface-portal-graph-v44';assert not O.exists();O.mkdir();C=G/'construction-contract.json';c=json.loads(C.read_bytes());E=G/'edge-search-evidence-v40.json';e=json.loads(E.read_bytes());ce=json.loads((G/'continued-search-evidence-v42.json').read_bytes());assert ce['completed']and not ce['route_found'];assert digest_file(G/'continued-progress-v42.npz')==ce['progress_sha256'];source=ROOT/c['source_catalogue'];assert digest_file(source)==c['source_catalogue_sha256'];atlas=json.loads(source.read_bytes());parts={p['id']:p for p in atlas['parts']};reader=AtlasReader();frame='model:local-ta98:assembled-metres';threshold=c['radius_m']+c['margin_m'];sg=reader.read_part(atlas,parts['TA98TONGUE:male:styloglossus:left'],frame=frame,catalogue_path=source,catalogue_sha=c['source_catalogue_sha256']);sp=reader.read_part(atlas,parts['ZA:Stylopharyngeus muscle.l'],frame=frame,catalogue_path=source,catalogue_sha=c['source_catalogue_sha256']);source_index={r['part_id']:i for i,r in enumerate(e['source_bindings'])};assert sg.geometry_sha==e['source_bindings'][source_index['TA98TONGUE:male:styloglossus:left']]['geometry_sha256'];assert sp.geometry_sha==e['source_bindings'][source_index['ZA:Stylopharyngeus muscle.l']]['geometry_sha256']
    contract=dict(schema='native-palatine-surface-portal-graph-contract/v44',source_catalogue_sha256=c['source_catalogue_sha256'],prior_terminal_graph_sha256=digest_file(G/'continued-search-evidence-v42.json'),source_count=179,radius_m=c['radius_m'],margin_m=c['margin_m'],SG_geometry_sha256=sg.geometry_sha,SP_geometry_sha256=sp.geometry_sha,
      method='Sample actual SG triangle barycentric surface at <=0.6 mm intervals, project to actual SP triangles; retain opposing-surface pairs with finite 1.3–8 mm separation. Sample 0.35/0.5/0.65 along these actual surface witness pairs; accept only full-radius points with real surface distance and valid closed-source occupancy. Build local 1.5 mm-neighbor graph, certify every bridge by adaptive exact source-triangle distance bounds. Reuse saved source-reachable native grid nodes and actual179-source field. This is not arbitrary dense local delta-search around the rejected tube.',
      Y_interval_m=[1.635,1.668],opposing_interval_separation_m=[.0013,.008],cross_section_fractions=[.35,.5,.65],maximum_node_neighbor_distance_m=.0015,source_grid_connector_distance_m=.0015,
      unchanged_source_station=8,unchanged_SG_SP_band=c['unchanged_SG_SP_band_m'],unchanged_recipient=c['unchanged_recipient'],fixed_other_routes=c['unchanged_other_route_hashes'],no_source_or_calibre_change=True,no_field_regeneration=True,no_C2_emitted_geometry_pack_install=True,internal_deadline_seconds=49)
    atomic_json(O/'construction-contract.json',contract)
    # Surface sampling is bounded by actual native triangles, independent of
    # the failed route or rejected guide's local displacements.
    samples=[];sample_faces=[]
    for fi,t in enumerate(sg.triangles):
        if t[:,1].max()<1.635 or t[:,1].min()>1.668:continue
        subdivisions=int(np.ceil(max(np.linalg.norm(t[1]-t[0]),np.linalg.norm(t[2]-t[0]),np.linalg.norm(t[2]-t[1]))/.0006));subdivisions=max(1,subdivisions)
        for i in range(subdivisions+1):
            for j in range(subdivisions+1-i):samples.append(t[0]+(t[1]-t[0])*(i/subdivisions)+(t[2]-t[0])*(j/subdivisions));sample_faces.append(fi)
    samples=np.array(samples);_,unique=np.unique(samples,axis=0,return_index=True);samples=samples[unique];sample_faces=np.array(sample_faces)[unique];d2,spfaces,nearest=igl.point_mesh_squared_distance(samples,sp.vertices,sp.faces);dist=np.sqrt(d2);selected=np.flatnonzero((dist>=.0013)&(dist<=.008)&(samples[:,1]>=1.635)&(samples[:,1]<=1.668));a=samples[selected];b=nearest[selected];points=np.concatenate([a+f*(b-a)for f in [.35,.5,.65]]);sg_faces=np.tile(sample_faces[selected],3);sp_faces=np.tile(spfaces[selected],3);fraction=np.repeat([.35,.5,.65],len(a));points,unique=np.unique(points,axis=0,return_index=True);sg_faces=sg_faces[unique];sp_faces=sp_faces[unique];fraction=fraction[unique]
    # Native opposing face vectors must place the point between both closest
    # surfaces rather than simply midway between unrelated geometry.
    da,fa,pa=igl.point_mesh_squared_distance(points,sg.vertices,sg.faces);db,fb,pb=igl.point_mesh_squared_distance(points,sp.vertices,sp.faces);na=(points-pa)/np.maximum(np.sqrt(da)[:,None],1e-30);nb=(points-pb)/np.maximum(np.sqrt(db)[:,None],1e-30);opposed=np.einsum('ij,ij->i',na,nb)<-.3;points=points[opposed];sg_faces=fa[opposed];sp_faces=fb[opposed];fraction=fraction[opposed]
    bank=np.load(G/'literal-surface-query-bank-v40.npz');assert digest_file(G/'literal-surface-query-bank-v40.npz')==e['query_bank_sha256'];tree=igl.AABB();tree.init(bank['vertices'],bank['faces']);d2,fi,near=tree.squared_distance(bank['vertices'],bank['faces'],points);distance=np.sqrt(d2);safe=distance>=threshold;classifiers=[]
    spec=importlib.util.spec_from_file_location('v44_occupancy',ROOT/'server/check-facial-cervical-feasibility.py');h=importlib.util.module_from_spec(spec);spec.loader.exec_module(h);B=D.parent;plan=json.loads((B/'male-left-facial-cervical-routes-cervical-c2-v17/routes.json').read_bytes());occupied=[];complete=True
    for row in e['source_bindings']:
        if time.monotonic()-began>45:complete=False;break
        pid=row['part_id']
        if row['classifier']not in['literal_closed_source','temporary_native_conduit_end_caps_never_rendered']:continue
        if pid in parts:m=reader.read_part(atlas,parts[pid],frame=frame,catalogue_path=source,catalogue_sha=c['source_catalogue_sha256']);allow=parts[pid].get('system')in{'arterial','venous'}
        elif pid=='CANDIDATE:glandular-0:passed-fixed-envelope':m=Mesh.load(ROOT/plan['extra_fixed_meshes'][0]['mesh']);allow=True
        else:
            name=pid.split(':')[1];z=np.load(D/(name+'-fields.npz'));m=Mesh(z['fields'][0],z['faces'],frame);allow=True
        assert m.geometry_sha==row['geometry_sha256'];cl,mode=h.occupancy_classifier(m,allow);assert cl is not None and mode==row['classifier'];ids=np.flatnonzero(safe&np.all(points>=cl.bounds[0],axis=1)&np.all(points<=cl.bounds[1],axis=1))
        if len(ids):inside=abs(igl.winding_number(cl.vertices,cl.faces,points[ids]))>.99;safe[ids[inside]]=False;occupied.append(dict(part_id=pid,tested=len(ids),inside_count=int(inside.sum())))
    atomic_npz(O/'source-surface-candidates.npz',points=points,safe=safe,source_sg_face=sg_faces,source_sp_face=sp_faces,surface_fraction=fraction,nearest_source_distance_m=distance,nearest_combined_face=fi)
    if not complete:atomic_json(O/'evidence.json',dict(completed=False,reason='node occupancy internal deadline',candidate_nodes=len(points),safe_nodes=int(safe.sum()),no_geometry_or_route=True,elapsed_seconds=time.monotonic()-began));return
    # Parent and PP are excluded from composite distance only for the original
    # prospective interfaces. This local interval is remote from those masks;
    # certify their surfaces explicitly for every candidate and edge.
    fixed=[]
    for pid,name in [('ZA:Facial artery.l','facial-left.npz'),('ZA:Palatopharyngeus muscle.l','palatopharyngeus-left.npz')]:
        m=Mesh.load(B/'male-left-facial-cervical-routes-anatomical-v5'/name);t=igl.AABB();t.init(m.vertices,m.faces);d2,_,_=t.squared_distance(m.vertices,m.faces,points);dd=np.sqrt(d2);safe&=dd>=threshold;distance=np.minimum(distance,dd);fixed.append((pid,t,m.vertices,m.faces))
    candidate_ids=np.flatnonzero(safe);new=points[candidate_ids];new_clear=distance[candidate_ids];grid=np.load(I/'clearance-domain.npz')['points'];progress=np.load(G/'continued-progress-v42.npz');reached=np.flatnonzero(np.isfinite(progress['distances'][0]));endband=np.array(c['unchanged_SG_SP_band_m']);target=np.flatnonzero(cKDTree(endband).query(new)[0]<=.002);newtree=cKDTree(new);gridtree=cKDTree(grid[reached]);connections=gridtree.query_ball_point(new,.0015);queue=[];cost=np.full(len(new),np.inf);prev=np.full(len(new),-1,dtype=np.int64);seeds=[];certs=[];edgecache={};stats=dict(calls=0,points=0,failed=0,unresolved=0)
    trees=[('composite',tree,bank['vertices'],bank['faces'])]+fixed
    def certify(a,b,lower):
        stats['calls']+=1;length=float(np.linalg.norm(b-a))
        if lower-length/2>=threshold:return True,dict(method='endpoint-Lipschitz',lower_m=lower-length/2)
        proof=[]
        for pid,t,v,f in trees:
            stack=[(0.,1.,0)]
            while stack:
                lo,hi,depth=stack.pop();mid=(lo+hi)/2;p=a+mid*(b-a);d2,fi,near=t.squared_distance(v,f,p[None]);dd=float(np.sqrt(d2[0]));lb=dd-(hi-lo)*length/2-3e-13;stats['points']+=1
                if lb>=threshold:proof.append(dict(surface=pid,interval=[lo,hi],distance_m=dd,lower_m=lb,face=int(fi[0])));continue
                if dd<threshold:stats['failed']+=1;return False,dict(method='actual-clearance-deficit',surface=pid,parameter=mid,distance_m=dd,face=int(fi[0]),source_point_m=near[0].tolist())
                if depth==12:stats['unresolved']+=1;return False,dict(method='depth-limit',surface=pid,parameter=mid,distance_m=dd,lower_m=lb,face=int(fi[0]))
                stack.extend([(lo,mid,depth+1),(mid,hi,depth+1)])
        return True,dict(method='adaptive-continuous-Lipschitz',proof=proof)
    field=json.loads((G/'field-extension-evidence.json').read_bytes())['current_field'];lower_grid=np.load(ROOT/field['path'])['clearance_lower_bound_m'];seed_attempts=0
    for i,nearby in enumerate(connections):
        if time.monotonic()-began>49:complete=False;break
        if not nearby:continue
        ds=np.linalg.norm(grid[reached[nearby]]-new[i],axis=1);order=np.argsort(ds)
        for j in order[:5]:
            gi=int(reached[nearby[j]]);ok,proof=certify(grid[gi],new[i],float(min(lower_grid[gi],new_clear[i])));seed_attempts+=1;certs.append(dict(grid_node=gi,local_node=i,passed=ok,proof=proof))
            if ok:
                nc=float(progress['distances'][0,gi]+ds[j]);cost[i]=nc;prev[i]=-gi-2;heapq.heappush(queue,(nc,i));seeds.append(i);break
    goal=set(map(int,target));found=None;pops=0
    while queue and complete:
        cc,i=heapq.heappop(queue)
        if cc!=cost[i]:continue
        pops+=1
        if time.monotonic()-began>49:complete=False;break
        if i in goal:found=i;break
        for j in newtree.query_ball_point(new[i],.0015):
            if j==i:continue
            pair=(min(i,j),max(i,j))
            if pair not in edgecache:
                ok,proof=certify(new[pair[0]],new[pair[1]],float(min(new_clear[list(pair)])));edgecache[pair]=ok;certs.append(dict(local_nodes=list(pair),passed=ok,proof=proof))
            if not edgecache[pair]:continue
            nc=cc+float(np.linalg.norm(new[i]-new[j]))
            if nc<cost[j]:cost[j]=nc;prev[j]=i;heapq.heappush(queue,(nc,j))
    atomic_npz(O/'local-graph-progress.npz',candidate_ids=candidate_ids,points=new,clearance_m=new_clear,distance=cost,previous=prev,queue=np.array(queue).reshape(-1,2),edges=np.array(list(edgecache),dtype=np.int64).reshape(-1,2),edge_pass=np.array(list(edgecache.values()),bool))
    route=None
    if found is not None:
        route=[];i=found
        while i>=0:route.append(i);i=int(prev[i])
        route=route[::-1];gi=-i-2;atomic_npz(O/'certified-source-to-band-connection.npz',grid_source_node=np.array(gi),local_node_ids=np.array(route),points=np.vstack([grid[gi],new[route]]),frame=np.array(frame))
    result=dict(schema='native-source-surface-adaptive-portal-connection/v44',source_catalogue_sha256=c['source_catalogue_sha256'],contract_sha256=digest_file(O/'construction-contract.json'),completed=complete,source_surface_samples=len(samples),opposing_candidate_nodes=len(points),actual179_source_clear_nodes=len(new),closed_occupancy_queries=occupied,source_grid_seed_attempts=seed_attempts,certified_source_grid_seeds=len(seeds),SG_SP_band_target_nodes=len(goal),visited_local_nodes=pops,local_edge_count=len(edgecache),statistics=stats,source_to_SG_SP_band_connection_found=found is not None,connection_sha256=digest_file(O/'certified-source-to-band-connection.npz')if found is not None else None,certificates=certs,progress_sha256=digest_file(O/'local-graph-progress.npz'),outcome='Source-surface diagnostic connection found; full path/C2/wall acceptance absent'if found is not None else 'No local source-surface connection certified in this finite sample; not anatomical impossibility',no_emitted_geometry_pack_install=True,elapsed_seconds=time.monotonic()-began)
    atomic_json(O/'evidence.json',result);print(json.dumps({k:result[k]for k in ['completed','source_surface_samples','opposing_candidate_nodes','actual179_source_clear_nodes','certified_source_grid_seeds','SG_SP_band_target_nodes','visited_local_nodes','statistics','source_to_SG_SP_band_connection_found','elapsed_seconds']}))

if __name__=='__main__':main()
