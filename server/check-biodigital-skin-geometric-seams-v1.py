"""Complete bounded segment/triangle tolerance intervals for residual seams.

No samples stand in for interval coverage; no source geometry or caps change.
All rendered-part triangle AABBs are searched, including large triangles that
a centroid/vertex shortlist could miss. The unused metric-hair picker carrier
is excluded from rendered-surface coverage.
"""
from pathlib import Path
from collections import defaultdict,Counter
import hashlib,json,sys,time
import numpy as np
start=time.monotonic();R=Path(__file__).resolve().parents[1];S=R/'.local-models/geometry-workbench/biodigital-male-circulation-v1/global-envelope-v1/fixed-skin';T=S/'envelope-topology-v1';O=T/'geometric-seam-coverage-v1'
if O.exists():raise SystemExit('Refusing to overwrite geometric-seam-coverage-v1')
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest();bindings={name:sha(S/name)for name in ['skin.npz','active-skin-layer.npz','policy.json','ownership.json']};skin=np.load(S/'skin.npz');p=skin['vertices'];f=skin['faces'];fo=skin['face_owner_part_index'];fs=skin['face_owner_source_face_id'];sv=skin['vertex_owner_source_vertex_id'];full=np.load(S/'active-skin-layer.npz');tp=full['vertices'];tf=full['faces'];towner=full['face_owner_part_index'];tsource=full['face_owner_source_face_id'];owners=json.loads((S/'ownership.json').read_bytes())['partOwnership'];by={r['partIndex']:r for r in owners};primary=np.array([i for i,r in by.items()if r['primaryBodySkin']],int);excluded=np.array([i for i,r in by.items()if' hairField' in r['sourceRecord']],int)
# Actual renderer replaces a hairField carrier with shaft instances.
rendered=np.array([i for i,r in by.items()if' hairField'.strip()not in r['sourceRecord']],int);ids=np.flatnonzero(np.isin(towner,rendered));tri=tp[tf[ids]];triowner=towner[ids];tri_low=tri.min(1);tri_high=tri.max(1);tri_primary=np.isin(triowner,primary);face_to_full={(int(a),int(b)):i for i,(a,b)in enumerate(zip(towner,tsource))};source_full=np.array([face_to_full[(int(a),int(b))]for a,b in zip(fo,fs)],int)
old=np.load(T/'proxy-and-source-evidence.npz');old_report=json.loads((T/'report.json').read_bytes());be=old['tol_1um_boundary_edges'];bface=old['tol_1um_boundary_original_skin_face_ids'];q=old['tol_1um_proxy_faces'];assert len(be)==1042;component_of={}
for c in old_report['toleranceLevels'][-1]['topology']['remainingBoundaryComponents']:
 for vertex in old[c['proxyVertexArray']]:component_of[int(vertex)]=c['componentIndex']
tolerance=1e-6;parameter_margin=2e-12;max_tests=200000;tested=0;saved={};witness=[];edge_rows=[];endpoint_rows=[]
def merge(intervals):
 if not intervals:return[]
 items=sorted((max(0.,float(a)),min(1.,float(b)))for a,b in intervals if b>=a and b>=0 and a<=1);result=[]
 for lo,hi in items:
  if result and lo<=result[-1][1]+parameter_margin:result[-1][1]=max(result[-1][1],hi)
  else:result.append([lo,hi])
 return result
def clip_linear(interval,a,b,lower=0.):
 if interval is None:return None
 lo,hi=interval
 if abs(b)<1e-30:return(lo,hi)if a>=lower-1e-14 else None
 root=(lower-a)/b
 if b>0:lo=max(lo,root)
 else:hi=min(hi,root)
 return(lo,hi)if lo<=hi+parameter_margin else None
def quadratic_ball(w,v,radius):
 a=float(v@v);b=float(2*w@v);c=float(w@w-radius*radius)
 if a<1e-32:return(0.,1.)if c<=0 else None
 disc=b*b-4*a*c;error=64*np.finfo(float).eps*(b*b+abs(4*a*c))
 if disc<-error:return None
 disc=max(0.,disc);s=np.sqrt(disc);lo=(-b-s)/(2*a);hi=(-b+s)/(2*a);return(max(0.,lo),min(1.,hi))if hi>=0 and lo<=1 else None
def triangle_intervals(A,B,t):
 # Local origin prevents cancellation of metre-scale body coordinates.
 x=t-A;D=B-A;e0=x[1]-x[0];e1=x[2]-x[0];n=np.cross(e0,e1);n2=float(n@n);intervals=[]
 if n2>0:
  nl=np.sqrt(n2);w=-x[0];l1=float(np.cross(w,e1)@n/n2);dl1=float(np.cross(D,e1)@n/n2);l2=float(np.cross(e0,w)@n/n2);dl2=float(np.cross(e0,D)@n/n2);iv=(0.,1.)
  for a,b in [(l1,dl1),(l2,dl2),(1-l1-l2,-dl1-dl2)]:iv=clip_linear(iv,a,b)
  d0=float(w@n/nl);dd=float(D@n/nl);iv=clip_linear(iv,tolerance-d0,-dd);iv=clip_linear(iv,tolerance+d0,dd)
  if iv:intervals.append(iv)
 for c,d in [(x[0],x[1]),(x[1],x[2]),(x[2],x[0])]:
  E=d-c;ll=float(E@E)
  if ll>0:
   s0=float((-c)@E/ll);sd=float(D@E/ll);w=-c-s0*E;v=D-sd*E;iv=quadratic_ball(w,v,tolerance);iv=clip_linear(iv,s0,sd);iv=clip_linear(iv,1-s0,-sd)
   if iv:intervals.append(iv)
  for vertex in [c,d]:
   iv=quadratic_ball(-vertex,D,tolerance)
   if iv:intervals.append(iv)
 return merge(intervals)
def stats(intervals,length):
 intervals=merge(intervals);total=sum(b-a for a,b in intervals);last=0.;gaps=[]
 for a,b in intervals:
  if a>last:gaps.append((last,a))
  last=max(last,b)
 if last<1:gaps.append((last,1.))
 maximum=max((b-a for a,b in gaps),default=0.);return dict(intervals=intervals,coveredFraction=float(total),coveredLength_m=float(total*length),uncoveredLength_m=float((1-total)*length),maximumUncoveredGap_m=float(maximum*length),fullCoverage=maximum<=1e-10)
def opposite_projection_interval(A,B,C,t):
 D=B-A;length=float(np.linalg.norm(D));u=D/length;inward=C-A-u*float((C-A)@u);il=float(np.linalg.norm(inward))
 if il==0:return None,0.
 inward/=il;points=[point.copy()for point in t];values=[float((point-A)@inward+2*tolerance)for point in points];clipped=[]
 for i,(point,value)in enumerate(zip(points,values)):
  prev=points[i-1];pv=values[i-1]
  if(value<=0)!=(pv<=0):clipped.append(prev+pv/(pv-value)*(point-prev))
  if value<=0:clipped.append(point)
 if not clipped:return None,float(min(values)-2*tolerance)
 pos=np.array([(point-A)@u/length for point in clipped]);return(float(pos.min()),float(pos.max())),float(min(values)-2*tolerance)
for edge_index,(edge,faceid)in enumerate(zip(be,bface)):
 faceid=int(faceid);corners=f[faceid];local=q[faceid];matches=[]
 for a,b in [(0,1),(1,2),(2,0)]:
  if np.array_equal(np.sort(local[[a,b]]),edge):matches.append((a,b))
 assert len(matches)==1;ia,ib=matches[0];ic=3-ia-ib;A=p[corners[ia]];B=p[corners[ib]];C=p[corners[ic]];D=B-A;length=float(np.linalg.norm(D));normal=np.cross(B-A,C-A);normal/=np.linalg.norm(normal);lo=np.minimum(A,B)-tolerance;hi=np.maximum(A,B)+tolerance;source_owner=int(fo[faceid]);source_global=int(source_full[faceid]);candidates=np.flatnonzero(np.all(tri_high>=lo,axis=1)&np.all(tri_low<=hi,axis=1)&(ids!=source_global));other=candidates[triowner[candidates]!=source_owner];same=candidates[triowner[candidates]==source_owner];complete=True
 if tested+len(candidates)>max_tests:complete=False;candidates=candidates[:max(0,max_tests-tested)]
 tested+=len(candidates);intervals=defaultdict(list);witness_count=0;target_counts=Counter()
 for target_index in candidates:
  t=tri[target_index];ivs=triangle_intervals(A,B,t)
  if not ivs:continue
  target_face=int(ids[target_index]);target_owner=int(triowner[target_index]);other_owner=target_owner!=source_owner;is_primary=bool(tri_primary[target_index]);tn=np.cross(t[1]-t[0],t[2]-t[0]);tnorm=float(np.linalg.norm(tn));cos=float(tn@normal/tnorm)if tnorm else None;opposite,minside=opposite_projection_interval(A,B,C,t);support=[]
  if opposite and cos is not None and cos>.5:
   support=merge([(max(a,opposite[0]),min(b,opposite[1]))for a,b in ivs if min(b,opposite[1])>=max(a,opposite[0])])
  scope='other'if other_owner else'same_owner_nonincident';intervals[scope+'_rendered'].extend(ivs)
  if is_primary:intervals[scope+'_primary'].extend(ivs)
  if other_owner:
   intervals['opposite_supported_other_rendered'].extend(support)
   if is_primary:intervals['opposite_supported_other_primary'].extend(support)
  for a,b in ivs:witness.append([edge_index,faceid,int(fs[faceid]),source_owner,target_face,int(tsource[target_face]),target_owner,a,b,int(other_owner),int(is_primary),cos if cos is not None else np.nan,minside]);witness_count+=1
  target_counts[by[target_owner]['id']]+=1
 row=dict(edgeIndex=edge_index,componentIndex=component_of[int(edge[0])],skinFaceId=faceid,sourcePartIndex=source_owner,sourcePartId=by[source_owner]['id'],originalSourceFaceId=int(fs[faceid]),originalSkinVertexIds=[int(corners[ia]),int(corners[ib])],originalSourceVertexIds=[int(sv[corners[ia]]),int(sv[corners[ib]])],rawEndpoints_m=[A.tolist(),B.tolist()],length_m=length,allTriangleAabbEnumerationComplete=complete,otherOwnerAabbCandidates=len(other),sameOwnerNonincidentAabbCandidates=len(same),triangleCandidatesTested=len(candidates),witnessIntervals=witness_count,coveringSourceParts=dict(target_counts),coverage={key:stats(intervals[key],length)for key in ['other_primary','other_rendered','opposite_supported_other_primary','opposite_supported_other_rendered','same_owner_nonincident_primary','same_owner_nonincident_rendered']});row['coverage']['all_nonincident_rendered']=stats(intervals['other_rendered']+intervals['same_owner_nonincident_rendered'],length);edge_rows.append(row);endpoint_rows.append([A,B])
 if not complete:break
component_rows=[]
for ci in sorted({r['componentIndex']for r in edge_rows}):
 rows=[r for r in edge_rows if r['componentIndex']==ci];original=next(c for c in old_report['toleranceLevels'][-1]['topology']['remainingBoundaryComponents']if c['componentIndex']==ci);scopes={}
 for key in edge_rows[0]['coverage']:
  scopes[key]=dict(fullEdges=sum(r['coverage'][key]['fullCoverage']for r in rows),testedEdges=len(rows),coveredLength_m=sum(r['coverage'][key]['coveredLength_m']for r in rows),uncoveredLength_m=sum(r['coverage'][key]['uncoveredLength_m']for r in rows),maximumUncoveredGap_m=max(r['coverage'][key]['maximumUncoveredGap_m']for r in rows))
 component_rows.append(dict(componentIndex=ci,sourceParts=original['sourceParts'],bounds_m=original['bounds_m'],simpleBoundaryCycle=original['simpleClosedBoundaryCycle'],coverage=scopes))
summary={}
for key in edge_rows[0]['coverage']:
 summary[key]=dict(fullEdges=sum(r['coverage'][key]['fullCoverage']for r in edge_rows),testedEdges=len(edge_rows),coveredLength_m=sum(r['coverage'][key]['coveredLength_m']for r in edge_rows),uncoveredLength_m=sum(r['coverage'][key]['uncoveredLength_m']for r in edge_rows),maximumUncoveredGap_m=max(r['coverage'][key]['maximumUncoveredGap_m']for r in edge_rows))
O.mkdir();np.savez_compressed(O/'interval-witnesses.npz',rows=np.array(witness,float).reshape(-1,13),raw_source_edge_endpoints_m=np.array(endpoint_rows),original_primary_boundary_edges=be,original_primary_boundary_source_faces=bface)
report=dict(schema='actual-source-geometric-seam-interval-coverage/v1',inputBindings=bindings,geometricTolerance_m=tolerance,intervalAlgorithm='Union of exact projected triangle-face/plane prism, finite edge-cylinder and vertex-sphere distance intervals for the complete line segment, in local coordinates; Float64 numeric predicates.',samplesUsedInsteadOfIntervals=False,candidateEnumeration='All actual rendered Skin-layer triangle AABBs, expanded by tolerance; excludes only unused hairField carrier; original owning face excluded. Other-owner and same-owner nonincident coverage reported separately.',targetTriangleCount=len(tri),maximumTriangleTests=max_tests,triangleCandidatesTested=tested,all1042EdgesTested=len(edge_rows)==1042,allAabbCandidatesTested=all(r['allTriangleAabbEnumerationComplete']for r in edge_rows),witnessColumns=['edgeIndex','skinFaceId','originalSourceFaceId','sourcePartIndex','allSkinTargetFace','targetOriginalSourceFace','targetSourcePartIndex','t0','t1','differentOwner','primaryTarget','faceNormalCosine','minimumOtherTriangleInteriorSide_m'],oppositeSupportQualification='Additional conservative diagnostic: target extends 2 micrometres to the opposite side of the owner face interior in projection and face normals agree within60degrees. Ordinary coverage alone may be an overlap; opposite support alone is not a complete closed-volume proof.',summary=summary,components=component_rows,edges=edge_rows,sourcePairMappingsRetained=True,absenceIsNotBasedOnShortlistMiss=True,unusedMetricHairCarrierNotCountedAsRenderedSkin=True,anatomicalCapsIntroduced=False,sourceGeometryChanged=False,globalInsideClassifierCertified=False,qualification='Geometric boundary coverage only. It does not by itself establish a closed oriented non-self-intersecting exterior union or account for material fragment-discard apertures. Ordinary hair-part geometry is identified separately from primary targets; hair shaft tips are not body clearance.',elapsedSeconds=time.monotonic()-start);(O/'report.json').write_text(json.dumps(report,indent=2)+'\n');assert all(sha(S/name)==digest for name,digest in bindings.items());print(json.dumps(dict(summary=summary,selectedComponents=[dict(componentIndex=c['componentIndex'],coverage=c['coverage'])for c in component_rows if c['componentIndex']in[41,5,79,17,54]],tested=tested,allEdges=report['all1042EdgesTested'],allCandidates=report['allAabbCandidatesTested'],elapsedSeconds=report['elapsedSeconds'])))
