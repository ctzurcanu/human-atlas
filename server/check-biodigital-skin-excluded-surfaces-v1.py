"""Compare existing nail/hair carrier coverage, without adding any geometry."""
from pathlib import Path
from collections import defaultdict,Counter
import hashlib,json,sys,time
import numpy as np
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components
start=time.monotonic();R=Path(__file__).resolve().parents[1];S=R/'.local-models/geometry-workbench/biodigital-male-circulation-v1/global-envelope-v1/fixed-skin';O=S/'envelope-topology-v1/excluded-surface-check-v1'
if O.exists():raise SystemExit('Refusing to overwrite excluded-surface-check-v1')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();bindings={n:sha(S/n)for n in ['skin.npz','active-skin-layer.npz','skin-appendages.npz','policy.json','ownership.json']};z=np.load(S/'active-skin-layer.npz');p=z['vertices'];f=z['faces'];fo=z['face_owner_part_index'];fs=z['face_owner_source_face_id'];vo=z['vertex_owner_part_index'];vs=z['vertex_owner_source_vertex_id'];owners=json.loads((S/'ownership.json').read_bytes())['partOwnership'];by={r['partIndex']:r for r in owners};primary=np.array([i for i,r in by.items()if r['primaryBodySkin']],int);nails=np.array([i for i,r in by.items()if any(g.lower()=='nail plate'for g in r['sourceRecord'].get('groups',[]))],int);metric=np.array([i for i,r in by.items()if r['id']=='ZA:Scalp hair:metric'],int);assert len(nails)==4 and len(metric)==1
old=np.load(S/'envelope-topology-v1/proxy-and-source-evidence.npz');bp=old['tol_1um_proxy_vertices_m'];be=old['tol_1um_boundary_edges'];bf=old['tol_1um_boundary_original_skin_face_ids'];skin=np.load(S/'skin.npz');skin_fo=skin['face_owner_part_index'];skin_fs=skin['face_owner_source_face_id'];full_lookup={(int(a),int(b)):i for i,(a,b)in enumerate(zip(fo,fs))};skin_face_to_full=np.array([full_lookup[(int(a),int(b))]for a,b in zip(skin_fo,skin_fs)],int);saved={};results=[];selections=[('primary_484',primary),('primary_plus_nails_488',np.r_[primary,nails]),('primary_plus_metric_carrier_485',np.r_[primary,metric]),('primary_nails_metric_carrier_489',np.r_[primary,nails,metric]),('all_skin_layer_493',np.array(list(by),int))]
def weld(v,tol=1e-6):
 pairs=cKDTree(v).query_pairs(tol,output_type='ndarray');d=np.linalg.norm(v[pairs[:,0]]-v[pairs[:,1]],axis=1);pairs=pairs[np.lexsort((pairs[:,1],pairs[:,0],d))];parent=np.arange(len(v));members={i:[i]for i in range(len(v))}
 def find(i):
  while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
  return int(i)
 for a,b in pairs:
  a=find(int(a));b=find(int(b))
  if a==b:continue
  if np.linalg.norm(v[members[a]][:,None]-v[members[b]][None],axis=2).max()>tol+1e-15:continue
  if b<a:a,b=b,a
  parent[b]=a;members[a]+=members.pop(b)
 rep=np.array([find(i)for i in range(len(v))]);unique,inv=np.unique(rep,return_inverse=True);return v[unique],inv
def edge_data(q):
 d=np.concatenate([q[:,[0,1]],q[:,[1,2]],q[:,[2,0]]]);e,inv,c=np.unique(np.sort(d,axis=1),axis=0,return_inverse=True,return_counts=True);order=np.argsort(inv,kind='stable');offset=np.r_[0,np.cumsum(c[:-1])];face=np.tile(np.arange(len(q)),3);sign=np.where(d[:,0]<d[:,1],1,-1);sums=np.bincount(inv,weights=sign,minlength=len(e));return d,e,inv,c,order,offset,face,sums
def boundary_components(v,e):
 if not len(e):return []
 graph=coo_matrix((np.ones(2*len(e),np.uint8),(np.r_[e[:,0],e[:,1]],np.r_[e[:,1],e[:,0]])),shape=(len(v),len(v))).tocsr();_,component=connected_components(graph,directed=False);degree=np.bincount(e.ravel(),minlength=len(v));result=[]
 for ci in np.unique(component[np.unique(e)]):
  ids=np.flatnonzero(component[e[:,0]]==ci);vertices=np.unique(e[ids]);result.append(dict(edges=len(ids),vertices=len(vertices),simpleCycle=bool(np.all(degree[vertices]==2)),bounds_m=[v[vertices].min(0).tolist(),v[vertices].max(0).tolist()]))
 return result
# Exact carrier/skin triangle equality determines whether the carrier adds area.
ev,ei=np.unique(p,axis=0,return_inverse=True);ef=ei[f];primary_keys=defaultdict(list)
for i in np.flatnonzero(np.isin(fo,primary)):primary_keys[tuple(sorted(map(int,ef[i])))].append(int(i))
carrier_pairs=[];carrier_faces=np.flatnonzero(np.isin(fo,metric));unmatched=[]
for i in carrier_faces:
 matches=primary_keys.get(tuple(sorted(map(int,ef[i]))),[])
 if not matches:unmatched.append(int(i))
 for j in matches:carrier_pairs.append([int(i),j,int(fo[i]),int(fo[j]),int(fs[i]),int(fs[j])])
saved['metric_carrier_exact_primary_face_pairs']=np.array(carrier_pairs,int).reshape(-1,6);saved['metric_carrier_faces_without_exact_primary_counterpart']=np.array(unmatched,int)
groups={'all_primary_residual_edges':np.arange(len(be))}
for name,ci in [('crown_cycle',41),('right_hand_network',5),('left_hand_network',79),('right_sole_network',17),('left_sole_network',54)]:
 ids=old['tol_1um_boundary_component_'+str(ci)+'_proxy_vertex_ids'];groups[name]=np.flatnonzero(np.isin(be[:,0],ids)&np.isin(be[:,1],ids))
for label,selected in selections:
 ids=np.flatnonzero(np.isin(fo,selected));used=np.unique(f[ids]);v,first,vi=np.unique(p[used],axis=0,return_index=True,return_inverse=True);global_to_exact=np.full(len(p),-1,int);global_to_exact[used]=vi;q=global_to_exact[f[ids]];pv,mapping=weld(v);q=mapping[q];d,e,inv,count,order,offset,face,sums=edge_data(q);boundary=e[count==1];components=boundary_components(pv,boundary);edge_lookup={tuple(k):i for i,k in enumerate(e)};point_lookup={tuple(point):int(mapping[i])for i,point in enumerate(v)};base_to_proxy=np.array([point_lookup[tuple(point)]for point in bp]);group_stats={};closure_pairs=[]
 for group,edge_ids in groups.items():
  c=Counter();extras=Counter()
  for oldeid in edge_ids:
   endpoints=base_to_proxy[be[oldeid]];key=tuple(sorted(map(int,endpoints)));edge_index=edge_lookup.get(key)
   if endpoints[0]==endpoints[1]:c['collapsedByTolerance']+=1;continue
   if edge_index is None:c['noSameProxyEdge']+=1;continue
   owners_idx=order[offset[edge_index]:offset[edge_index]+count[edge_index]];global_faces=ids[face[owners_idx]];is_extra=~np.isin(fo[global_faces],primary);extra_faces=global_faces[is_extra]
   if count[edge_index]==1:c['stillBoundary']+=1
   elif count[edge_index]==2 and sums[edge_index]==0:
    c['manifoldOppositeEdge']+=1
    if len(extra_faces):
     c['properClosureByExistingAppendage']+=1
     for af in extra_faces:
      extras[by[int(fo[af])]['id']]+=1;body_face=int(skin_face_to_full[bf[oldeid]]);closure_pairs.append([int(oldeid),int(bf[oldeid]),body_face,int(af),int(fo[body_face]),int(fo[af]),int(fs[body_face]),int(fs[af])])
   elif count[edge_index]==2:c['sameWindingOverlayNotClosure']+=1
   else:c['nonmanifoldOverlayNotClosure']+=1
  group_stats[group]=dict(originalPrimaryBoundaryEdges=len(edge_ids),counts=dict(c),appendageSourceParts=dict(extras))
 saved[label+'_proper_closure_source_pairs']=np.unique(np.array(closure_pairs,int).reshape(-1,8),axis=0);saved[label+'_vertices_m']=pv;saved[label+'_faces']=q;saved[label+'_original_all_skin_face_ids']=ids;saved[label+'_boundary_edges']=boundary
 results.append(dict(label=label,parts=len(selected),faces=len(q),proxyTolerance_m=1e-6,boundaryEdges=len(boundary),boundaryComponents=len(components),simpleCycles=sum(c['simpleCycle']for c in components),branchedNetworks=sum(not c['simpleCycle']for c in components),nonmanifoldEdges=int(np.sum(count>2)),inconsistentWindingEdges=int(np.sum((count==2)&(np.abs(sums)>0))),duplicateFaces=int(len(q)-len(np.unique(np.sort(q,axis=1),axis=0))),zeroAreaFaces=int(np.sum(np.linalg.norm(np.cross(pv[q[:,1]]-pv[q[:,0]],pv[q[:,2]]-pv[q[:,0]]),axis=1)==0)),sourceClosureComparisons=group_stats,proxyOnly=True,closedBodyCertified=False))
# Distances identify actual source counterparts; distance alone is not closure.
proximity=[]
for group,edge_ids in groups.items():
 if group=='all_primary_residual_edges':continue
 points=bp[np.unique(be[edge_ids])]
 for pi in np.r_[nails,metric]:
  ap=p[np.unique(f[np.flatnonzero(fo==pi)])];dist,nearest=cKDTree(ap).query(points);key=group+'_part_'+str(pi);saved[key+'_primary_boundary_points_m']=points;saved[key+'_appendage_points_m']=ap[nearest];saved[key+'_source_part_index']=np.full(len(points),int(pi));saved[key+'_distance_m']=dist;proximity.append(dict(group=group,sourcePartIndex=int(pi),sourcePartId=by[int(pi)]['id'],boundaryVertices=len(points),minimumVertexPairDistance_m=float(dist.min()),medianVertexPairDistance_m=float(np.median(dist)),maximumVertexPairDistance_m=float(dist.max()),matchesWithin0p1Micrometre=int(np.sum(dist<=1e-7)),matchesWithin1Micrometre=int(np.sum(dist<=1e-6)),sourcePairArrayPrefix=key,distancesAloneNotClosureProof=True))
scene=R/'app/scene.tsx';hair=R/'app/metric-hair-field.ts';O.mkdir();np.savez_compressed(O/'source-and-closure-evidence.npz',**saved);report=dict(schema='existing-active-skin-appendage-closure-check/v1',inputBindings=bindings,rendererCodeHashes={str(scene.relative_to(R)):sha(scene),str(hair.relative_to(R)):sha(hair)},metricHairSemantics=dict(carrierUsedForPicking=True,carrierRenderedAsBodySurface=False,visibleGeometry='Closed instanced cylinder shafts, 1 mm long and 0.1 mm diameter; roots and normals from the hairField bank.',skinBetweenStrandsMustComeFromBodySkin=True,hairStrandTipsNotBodyClearance=True),metricCarrierCoverage=dict(sourcePart='ZA:Scalp hair:metric',carrierFaces=len(carrier_faces),carrierFacesWithExactPrimarySkinCounterpart=len(carrier_faces)-len(unmatched),carrierFacesWithoutExactPrimaryCounterpart=len(unmatched),pairArray='metric_carrier_exact_primary_face_pairs',pairColumns=['allSkinCarrierFace','allSkinPrimaryFace','carrierSourcePartIndex','primarySourcePartIndex','carrierOriginalSourceFace','primaryOriginalSourceFace']),unions=results,actualSourcePairDistances=proximity,sourcePartIds={str(i):r['id']for i,r in by.items()if r['skinLayerAccessory']},noSourceOrBodyChange=True,noAnatomicalCapsIntroduced=True,noHairTipsUsedForBodyClearance=True,elapsedSeconds=time.monotonic()-start);(O/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert all(sha(S/n)==digest for n,digest in bindings.items());print(json.dumps(dict(metricCarrierCoverage=report['metricCarrierCoverage'],unions=[{k:r[k]for k in ['label','boundaryEdges','simpleCycles','branchedNetworks','nonmanifoldEdges','inconsistentWindingEdges','duplicateFaces']}|{'crown':r['sourceClosureComparisons']['crown_cycle'],'rightHand':r['sourceClosureComparisons']['right_hand_network'],'rightSole':r['sourceClosureComparisons']['right_sole_network']}for r in results],elapsedSeconds=report['elapsedSeconds'])))
