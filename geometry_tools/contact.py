"""Actual near-surface patches, reproducible witnesses and bounded patch area."""
from pathlib import Path
from collections import defaultdict
import os,time
import numpy as np
from .common import CACHE,VERSION,identity,digest_file,atomic_npz,owned_path
from .distance import triangle_pairs,point_triangle
from .spatial import BVH,cached_bvh


def require_common_frame(first,second):
    if first.frame!=second.frame:
        raise ValueError('Different coordinate frames: register an explicit transform before comparing or fitting these meshes')


def pair_contacts(first,second,distance_m,*,first_bvh=None,second_bvh=None,selfcheck=False,max_pairs=250000,max_tests=2000000):
    require_common_frame(first,second)
    if not isinstance(max_pairs,int)or max_pairs<1:raise ValueError('Positive finite pair budget required')
    a=first_bvh or BVH(first);b=second_bvh or BVH(second);pieces=defaultdict(list);pairs=[];candidate_count=0;found=0;truncated=False
    if not isinstance(max_tests,int)or max_tests<1:raise ValueError('Positive finite triangle-test budget required')
    for candidate in a.candidate_batches(b,distance_m,selfcheck=selfcheck):
        if candidate_count+len(candidate)>max_tests:candidate=candidate[:max_tests-candidate_count];truncated=True
        if not len(candidate):break
        candidate_count+=len(candidate);r=triangle_pairs(a.triangles[candidate[:,0]],b.triangles[candidate[:,1]])
        take=r['distance_m']<=distance_m;ids=np.flatnonzero(take)
        if found+len(ids)>max_pairs:ids=ids[:max_pairs-found];truncated=True
        if len(ids):
            pairs.append(candidate[ids]);found+=len(ids)
            for key,value in r.items():pieces[key].append(value[ids])
        if truncated:break
    result={key:np.concatenate(values)for key,values in pieces.items()};result['pairs']=np.concatenate(pairs)if pairs else np.empty((0,2),dtype='<i8')
    if not pairs:
        empty=triangle_pairs(np.empty((0,3,3)),np.empty((0,3,3)));result.update(empty)
    summary=dict(distance_threshold_m=float(distance_m),candidate_triangle_pairs_tested=candidate_count,returned_triangle_pairs=found,truncated=truncated,
                 complete_within_pair_budget=not truncated,first_patch_faces=int(len(np.unique(result['pairs'][:,0]))),second_patch_faces=int(len(np.unique(result['pairs'][:,1]))),
                 crossing_face_pairs=int(np.sum(result['relation']=='crossing')),coplanar_overlap_face_pairs=int(np.sum(result['relation']=='coplanar-overlap')),
                 touching_face_pairs=int(np.sum(result['relation']=='touching')),degenerate_pair_count=int(np.sum(result['degenerate'])),
                 anatomical_attachment_inferred=False,anatomical_acceptance=False)
    for side,mesh in [('first',first),('second',second)]:
        col=0 if side=='first'else 1;ids=np.unique(result['pairs'][:,col]);t=mesh.triangles[ids]
        summary[side+'_candidate_patch_area_upper_bound_m2']=float(np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1).sum()/2)
    return result,summary


def self_contact_summary(mesh,*,bvh=None,root_pair=None,max_pairs=250000,max_tests=2000000,stop_on_defect=False,example_limit=20):
    """Same exact triangle tests, with streaming reduction of touching pairs.

    A closed wall has millions of legitimate adjacent contacts. A self audit
    needs their counts and face masks, not every barycentric witness in RAM.
    Explicit BVH root pairs also permit separately bounded saved stages.
    """
    if not isinstance(max_pairs,int)or max_pairs<1 or not isinstance(max_tests,int)or max_tests<1:raise ValueError('Positive finite pair/test budgets required')
    tree=bvh or BVH(mesh);first=np.zeros(len(mesh.faces),bool);second=first.copy();tested=found=defects=degenerate=0;counts=defaultdict(int);examples=[];truncated=False;stopped=False
    for candidate in tree.candidate_batches(tree,0.,selfcheck=True,root_pair=root_pair):
        if tested+len(candidate)>max_tests:candidate=candidate[:max_tests-tested];truncated=True
        if not len(candidate):break
        tested+=len(candidate);r=triangle_pairs(tree.triangles[candidate[:,0]],tree.triangles[candidate[:,1]]);ids=np.flatnonzero(r['distance_m']<=0.)
        if found+len(ids)>max_pairs:ids=ids[:max_pairs-found];truncated=True
        pairs=candidate[ids];relations=r['relation'][ids];found+=len(ids);degenerate+=int(r['degenerate'][ids].sum())
        first[pairs[:,0]]=True;second[pairs[:,1]]=True
        for name in ['crossing','coplanar-overlap','touching']:counts[name]+=int(np.sum(relations==name))
        wrong=np.isin(relations,['crossing','coplanar-overlap']);defects+=int(wrong.sum());examples.extend(pairs[wrong][:max(0,example_limit-len(examples))].tolist())
        if stop_on_defect and np.any(wrong):stopped=True;break
        if truncated:break
    summary=dict(distance_threshold_m=0.,candidate_triangle_pairs_tested=tested,returned_triangle_pairs=found,truncated=truncated,complete_within_pair_budget=not(truncated or stopped),
        first_patch_faces=int(first.sum()),second_patch_faces=int(second.sum()),crossing_face_pairs=counts['crossing'],coplanar_overlap_face_pairs=counts['coplanar-overlap'],touching_face_pairs=counts['touching'],degenerate_pair_count=degenerate,
        anatomical_attachment_inferred=False,anatomical_acceptance=False,defect_pairs=defects,examples=examples)
    area=np.linalg.norm(np.cross(mesh.triangles[:,1]-mesh.triangles[:,0],mesh.triangles[:,2]-mesh.triangles[:,0]),axis=1)/2
    summary['first_candidate_patch_area_upper_bound_m2']=float(area[first].sum());summary['second_candidate_patch_area_upper_bound_m2']=float(area[second].sum())
    return summary,dict(first=np.flatnonzero(first),second=np.flatnonzero(second))


def patch_components(mesh,face_ids):
    """Connected actual source-face patches, with no proximity-based welding."""
    ids=np.unique(np.asarray(face_ids,int));byvertex=defaultdict(list)
    # Attribute seams may duplicate controls: weld only literally equal points.
    v,iv=np.unique(mesh.vertices,axis=0,return_inverse=True);f=iv[mesh.faces]
    for face in ids:
        for vertex in f[face]:byvertex[int(vertex)].append(int(face))
    remaining=set(map(int,ids));groups=[]
    while remaining:
        stack=[min(remaining)];group=[]
        while stack:
            face=stack.pop()
            if face not in remaining:continue
            remaining.remove(face);group.append(face)
            for vertex in f[face]:stack.extend(byvertex[int(vertex)])
        groups.append(sorted(group))
    return groups


def patch_area_bounds(mesh,other_bvh,face_ids,distance_m,*,absolute_area_error_m2=1e-6,max_depth=9,max_subtriangles=20000):
    """Numerical lower/upper bounds from convexity and 1-Lipschitz distance.

    A face with minimum distance below the threshold is not automatically
    counted as a complete near-contact face. Remaining ambiguity is reported.
    """
    if distance_m<0 or absolute_area_error_m2<=0 or max_depth<0:raise ValueError('Invalid contact area controls')
    ids=np.unique(np.asarray(face_ids,int));triangles=mesh.triangles[ids]
    area=lambda t:float(np.linalg.norm(np.cross(t[1]-t[0],t[2]-t[0]))/2)
    total=sum(area(t)for t in triangles)
    if total==0:return dict(lower_m2=0.,upper_m2=0.,uncertainty_m2=0.,target_error_met=True,subtriangles_tested=0)
    pending=[(t,0,absolute_area_error_m2*area(t)/total)for t in triangles];lower=upper=0.;tested=0;unresolved=0
    while pending:
        if tested>=max_subtriangles:
            upper+=sum(area(t)for t,_,_ in pending);unresolved+=len(pending);break
        t,depth,budget=pending.pop();a=area(t);tested+=1
        center=t.mean(0);q=other_bvh.closest_point(center);radius=float(np.linalg.norm(t-center,axis=1).max())
        roundoff=512*np.finfo(float).eps*max(float(np.max(abs(t))),float(np.max(abs(other_bvh.low[0]))),float(np.max(abs(other_bvh.high[0]))),np.finfo(float).tiny)
        if q['distance_m']-radius-roundoff>distance_m:continue
        # Distance to one convex triangle is convex; if all three vertices
        # lie within its distance offset, the whole source subtriangle does.
        d,_,_=point_triangle(t,np.broadcast_to(other_bvh.triangles[q['face_id']],(3,3,3)))
        if np.sqrt(np.max(d))+roundoff<=distance_m:
            lower+=a;upper+=a;continue
        if q['distance_m']+radius+roundoff<=distance_m:
            lower+=a;upper+=a;continue
        if depth>=max_depth or a<=budget:
            upper+=a;unresolved+=1;continue
        x,y,z=(t[0]+t[1])*.5,(t[1]+t[2])*.5,(t[2]+t[0])*.5
        pending.extend((child,depth+1,budget/4)for child in [np.array([t[0],x,z]),np.array([x,t[1],y]),np.array([z,y,t[2]]),np.array([x,y,z])])
    return dict(lower_m2=lower,upper_m2=upper,uncertainty_m2=max(0.,upper-lower),target_error_m2=absolute_area_error_m2,
                target_error_met=upper-lower<=absolute_area_error_m2,subtriangles_tested=tested,unresolved_subtriangles=unresolved,work_budget_exhausted=bool(pending and tested>=max_subtriangles),
                bound_method='Float64 source subtriangle areas; convex-triangle offset and 1-Lipschitz bounds with scale-dependent roundoff margin',exact_arithmetic=False)


def contacts(db,first,second,distance_m,*,selfcheck=False,max_pairs=250000,area_error_m2=None,area_max_depth=8):
    require_common_frame(first,second)
    parameters=dict(distance_m=float(distance_m),selfcheck=selfcheck,max_pairs=max_pairs,area_error_m2=area_error_m2,area_max_depth=area_max_depth,
                    contact_implementation=digest_file(__file__),distance_implementation=digest_file(Path(__file__).with_name('distance.py')))
    key=identity(VERSION,'contact',first.geometry_sha,second.geometry_sha,first.frame,parameters);cached=db.cache_get(key)
    if cached:
        with np.load(cached['path'],allow_pickle=False)as z:arrays={k:z[k].copy()for k in z.files}
        return arrays,{**cached['summary'],'cache_hit':True,'cache_key':key,'artifact':cached['path']}
    began=time.monotonic();a,ac=cached_bvh(db,first);b,bc=cached_bvh(db,second)
    arrays,summary=pair_contacts(first,second,distance_m,first_bvh=a,second_bvh=b,selfcheck=selfcheck,max_pairs=max_pairs)
    summary['first_patch_components']=patch_components(first,arrays['pairs'][:,0]);summary['second_patch_components']=patch_components(second,arrays['pairs'][:,1])
    if area_error_m2 is not None and not summary['truncated']:
        summary['first_area_bounds']=patch_area_bounds(first,b,arrays['pairs'][:,0],distance_m,absolute_area_error_m2=area_error_m2,max_depth=area_max_depth)
        summary['second_area_bounds']=patch_area_bounds(second,a,arrays['pairs'][:,1],distance_m,absolute_area_error_m2=area_error_m2,max_depth=area_max_depth)
    summary.update(elapsed_seconds=time.monotonic()-began,bvh_cache=[ac,bc],frame=first.frame)
    path=CACHE/'contacts'/(key+'.npz')
    if path.exists():owned_path(path);os.replace(path,path.with_suffix('.corrupt-'+str(time.time_ns())+'.npz'))
    atomic_npz(path,**arrays);db.cache_put(key=key,geometry_sha=first.geometry_sha,kind='surface-contact',parameters=parameters,path=path,summary=summary)
    return arrays,{**summary,'cache_hit':False,'cache_key':key,'artifact':str(path)}
