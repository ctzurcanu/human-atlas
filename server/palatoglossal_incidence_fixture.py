"""Source-bound exact ancestry fixture; no production/core predicates change.

The Boolean only supplies topology/ancestry. Rational controls are reconstructed
from the original finite source triangles before comparing same-parent charts.
"""
from fractions import Fraction as Q
from itertools import combinations
import numpy as np
from geometry_tools.algorithms.exact_source_intersections import contains,edge_plane,_sub,_dot,_cross

def qvec(v):return tuple(Q.from_float(float(x))for x in v)

def rational_controls(raw,sources,maximum_replay_m=2e-10,selected_points=None):
    from geometry_tools.algorithms.finite_tissue_boolean import provenance
    v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts);oid,pid=provenance(raw)
    incident={}
    for face,original,parent in zip(f,oid,pid):
        for i in face:incident.setdefault(tuple(v[i]),set()).add((int(original),int(parent)))
    original_points={tuple(p):qvec(p)for source in sources.values()for p in source.vertices}
    result={};evidence={}
    for point,parents in incident.items():
        if selected_points is not None and point not in selected_points:continue
        if point in original_points:
            result[point]=original_points[point];continue
        tris=[tuple(qvec(p)for p in sources[original].vertices[sources[original].faces[parent]])for original,parent in sorted(parents)]
        candidates=set(q for tri in tris for q in tri)
        # Edges/planes suffice for this retained finite arrangement. A missing
        # valid triple-plane candidate is diagnosed separately, never guessed.
        for a,b in set(edge for tri in tris for edge in combinations(tri,2)):
            for triangle in tris:
                hit=edge_plane(a,b,triangle)
                if hit is not None:candidates.add(hit)
        valid=[]
        for candidate in candidates:
            if all(contains(candidate,t)for t in tris):
                d=float(np.linalg.norm(np.array([float(x)for x in candidate])-point))
                if d<=maximum_replay_m:valid.append((d,candidate))
        if not valid:raise ValueError('Original source ancestry is inconsistent for '+str((point,sorted(parents))))
        distance,q=min(valid,key=lambda x:x[0]);result[point]=q;evidence[point]=dict(parents=sorted(parents),replay_m=distance)
    return result,evidence

def plane_coordinates(triangle):
    n=_cross(_sub(triangle[1],triangle[0]),_sub(triangle[2],triangle[0]));drop=max(range(3),key=lambda i:abs(n[i]));return [i for i in range(3)if i!=drop]

def oriented_area(poly):
    return sum(a[0]*b[1]-a[1]*b[0]for a,b in zip(poly,poly[1:]+poly[:1]))/2

def clip_polygon(subject,cutter):
    """Exact rational convex polygon clipping in one original parent chart."""
    poly=list(subject);cut=list(cutter)
    if oriented_area(cut)<0:cut.reverse()
    for a,b in zip(cut,cut[1:]+cut[:1]):
        if not poly:break
        def side(p):return (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
        out=[]
        for p,q in zip(poly,poly[1:]+poly[:1]):
            sp,sq=side(p),side(q)
            if sp>=0:out.append(p)
            if (sp<0<sq)or(sq<0<sp):
                t=sp/(sp-sq);out.append(tuple(p[k]+t*(q[k]-p[k])for k in range(2)))
            elif sp<0 and sq==0:out.append(q)
        poly=list(dict.fromkeys(out))
    return poly

def common_parent_arrangement(triangles,parent_ids,original_triangle):
    if len(set(parent_ids))!=1:
        raise ValueError('Distinct original source planes must remain distinct')
    if not all(contains(p,original_triangle)for t in triangles for p in t):
        raise ValueError('A fragment has not been reconstructed on its literal original finite source face')
    axes=plane_coordinates(original_triangle)
    charts=[[tuple(p[i]for i in axes)for p in t]for t in triangles]
    common=charts[0]
    for chart in charts[1:]:common=clip_polygon(common,chart)
    return axes,common
