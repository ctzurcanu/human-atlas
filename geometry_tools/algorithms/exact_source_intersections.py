"""Replay Boolean cut controls on the literal binary source triangles.

The Boolean engine supplies incidence, not permission to distance-weld two
different features. Exact rational edge/plane intersections are evaluated once
and reused by every owner. Numerical replay distance is recorded separately.
"""
from fractions import Fraction as Q
from itertools import combinations
import numpy as np


def _sub(a,b):return tuple(x-y for x,y in zip(a,b))
def _dot(a,b):return sum(x*y for x,y in zip(a,b))
def _cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])


def edge_plane(first,last,triangle):
    a,b,c=triangle;n=_cross(_sub(b,a),_sub(c,a));d=_sub(last,first);den=_dot(n,d)
    if not den:return None
    t=_dot(n,_sub(a,first))/den
    if not 0<=t<=1:return None
    return tuple(x+t*y for x,y in zip(first,d))


def contains(point,triangle):
    a,b,c=triangle;u,v,w=_sub(b,a),_sub(c,a),_sub(point,a);n=_cross(u,v)
    if _dot(n,w):return False
    drop=max(range(3),key=lambda i:abs(n[i]));axes=[i for i in range(3)if i!=drop];i,j=axes
    det=u[i]*v[j]-u[j]*v[i]
    if not det:return False
    s=(w[i]*v[j]-w[j]*v[i])/det;t=(u[i]*w[j]-u[j]*w[i])/det
    return s>=0 and t>=0 and s+t<=1


def replay_boolean_vertices(vertices,faces,original_ids,parent_ids,
    first_id,first_positions,first_faces,second_id,second_positions,second_faces,
    *,maximum_replay_m=2e-10):
    v=np.asarray(vertices,float);f=np.asarray(faces,int)
    pp=[np.asarray(first_positions,float),np.asarray(second_positions,float)]
    ff=[np.asarray(first_faces,int),np.asarray(second_faces,int)]
    family={int(first_id):0,int(second_id):1};incident=[({}, {})for _ in v]
    exact={tuple(point)for p,fs in zip(pp,ff)for point in p[np.unique(fs)]}
    for tri,oid,pid in zip(f,original_ids,parent_ids):
        k=family[int(oid)]
        for index in tri:incident[int(index)][k][int(pid)]=None
    rationals={};triangles={}
    def rational(point):
        key=tuple(point)
        if key not in rationals:rationals[key]=tuple(Q.from_float(float(x))for x in key)
        return rationals[key]
    def triangle(k,pid):
        key=k,pid
        if key not in triangles:triangles[key]=tuple(rational(p)for p in pp[k][ff[k][pid]])
        return triangles[key]
    result=v.copy();distances=[];changed=0
    for index,point in enumerate(v):
        if tuple(point)in exact:continue
        parents=incident[index]
        source=[[triangle(k,pid)for pid in parent]for k,parent in enumerate(parents)]
        candidates=set()
        for group in source:
            if group:
                common=set(group[0])
                for tri in group[1:]:common.intersection_update(tri)
                if len(common)==1:candidates.update(common)
        if not all(source):
            candidates.update(q for group in source for tri in group for q in tri)
        else:
            for k in [0,1]:
                common=set(source[k][0])
                for tri in source[k][1:]:common.intersection_update(tri)
                edges=list(combinations(sorted(common),2))if len(common)==2 else [
                    (tri[i],tri[(i+1)%3])for tri in source[k]for i in range(3)]
                for a,b in edges:
                    for tri in source[1-k]:
                        hit=edge_plane(a,b,tri)
                        if hit is not None:candidates.add(hit)
        accepted=[]
        for q in candidates:
            if all(contains(q,tri)for group in source for tri in group):
                value=np.array([float(x)for x in q]);d=float(np.linalg.norm(value-point))
                if d<=maximum_replay_m:accepted.append((d,value))
        if not accepted:
            raise ValueError(f'No exact incident source intersection for output control {index}; '
                f'parents={[list(p)for p in parents]}, point={point.tolist()}')
        distance,value=min(accepted,key=lambda pair:pair[0])
        result[index]=value;distances.append(distance);changed+=1
    return result,dict(exact_source_intersection_controls=changed,
        maximum_source_replay_m=max(distances,default=0.),
        maximum_allowed_source_replay_m=maximum_replay_m,
        input_source_controls_literal=True,distance_welding=False)
