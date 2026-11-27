"""Producer-specific original finite-face and intermediate-event ancestry.

Original IDs/face IDs remain authoritative. A previously computed exact event
may survive after one creating face becomes interior to a later Boolean. Its
literal output-coordinate identity is retained with its rational provenance;
nearness or coplanarity never creates ancestry.
"""
from itertools import combinations
from fractions import Fraction as Q
import numpy as np
from collections import Counter
from geometry_tools.algorithms.exact_source_intersections import contains,edge_plane,_sub,_dot,_cross
from geometry_tools.algorithms.finite_tissue_boolean import provenance
from palatoglossal_incidence_fixture import qvec

def interior_steiner_control(point,star_triangles,original_parents,source_triangle):
    """Exact interior chart control; never a boundary/intersection exception.

    Only a complete closed, consistently oriented same-parent vertex star can
    take this producer-specific route. The rational control must be strictly
    interior to its literal finite parent and round to the identical binary
    point, so no geometric point or anatomical interface is moved or welded.
    """
    if len(set(original_parents))!=1:raise ValueError('Steiner control has multiple original parents')
    key=tuple(point);links=[]
    for triangle in star_triangles:
        corners=[tuple(p)for p in triangle];where=[i for i,p in enumerate(corners)if p==key]
        if len(where)!=1:raise ValueError('Steiner star has an invalid or degenerate center')
        i=where[0];links.append((corners[(i+1)%3],corners[(i+2)%3]))
    if len(links)<3 or len(set(links))!=len(links):raise ValueError('Steiner star is incomplete')
    outgoing=Counter(a for a,b in links);incoming=Counter(b for a,b in links)
    if set(outgoing)!=set(incoming)or any(n!=1 for n in outgoing.values())or any(n!=1 for n in incoming.values()):raise ValueError('Steiner star is not a closed oriented link')
    next_vertex=dict(links);start=links[0][0];walk=set();at=start
    while at not in walk:walk.add(at);at=next_vertex[at]
    if at!=start or len(walk)!=len(links):raise ValueError('Steiner star contains disconnected link cycles')
    t=tuple(qvec(p)for p in source_triangle);a,b,c=t;n=_cross(_sub(b,a),_sub(c,a));den=_dot(n,n)
    if not den:raise ValueError('Steiner parent is degenerate')
    qp=qvec(point);residual=_dot(n,_sub(qp,a));scale=residual/den;q=tuple(qp[i]-scale*n[i]for i in range(3))
    if tuple(float(x)for x in q)!=key:
        # Orthogonal projection is not the appropriate metric for an
        # anisotropic Float64 rounding box. Intersect that box with the exact
        # source plane analytically. This is one closed-form solve, never a
        # tolerance increase, geometric parameter trial or nearby-point weld.
        toward=[]
        for i in range(3):
            direction=-1 if residual*n[i]>0 else 1
            neighbor=np.nextafter(float(point[i]),-np.inf if direction<0 else np.inf)
            half=abs(Q.from_float(float(neighbor))-qp[i])/2
            toward.append(direction*half if n[i]else Q(0))
        support=-_dot(n,toward)if residual>0 else _dot(n,toward)
        if support<=0 or abs(residual)>support:raise ValueError('Steiner source plane misses the identical Float64 rounding cell')
        fraction=abs(residual)/support;q=tuple(qp[i]+fraction*toward[i]for i in range(3))
    if not contains(q,t):raise ValueError('Steiner control leaves its exact finite parent')
    if any(not any(_cross(_sub(t[(i+1)%3],t[i]),_sub(q,t[i])))for i in range(3)):raise ValueError('Steiner control is on an original parent boundary')
    if tuple(float(x)for x in q)!=key:raise ValueError('Steiner rational encoding would move a binary control')
    return q

def original_controls(raw,sources,inherited,*,stage,maximum_replay_m=2e-10):
    v=np.asarray(raw.vert_properties)[:,:3];f=np.asarray(raw.tri_verts);oid,pid=provenance(raw)
    incident={};stars={}
    for face,original,parent in zip(f,oid,pid):
        for i in face:
            point=tuple(v[i]);incident.setdefault(point,set()).add((int(original),int(parent)));stars.setdefault(point,[]).append(v[face])
    original_points={tuple(p):qvec(p)for source in sources.values()for p in source.vertices};result={};events=[]
    for point,parents in incident.items():
        triangles=[tuple(qvec(p)for p in sources[original].vertices[sources[original].faces[parent]])for original,parent in sorted(parents)]
        if point in inherited:
            candidates={inherited[point]['point']};method='literal-preserved-intermediate-event'
        elif point in original_points:
            candidates={original_points[point]};method='literal-source-control'
        else:
            candidates=set(p for triangle in triangles for p in triangle);method='original-finite-edge-plane-intersection'
            for a,b in set(edge for triangle in triangles for edge in combinations(triangle,2)):
                for triangle in triangles:
                    hit=edge_plane(a,b,triangle)
                    if hit is not None:candidates.add(hit)
        valid=[]
        for q in candidates:
            if all(contains(q,t)for t in triangles):
                d=float(np.linalg.norm(np.array([float(x)for x in q])-point))
                if d<=maximum_replay_m:valid.append((d,q))
        if not valid and len(parents)==1:
            parent_original,parent=next(iter(parents))
            try:q=interior_steiner_control(point,stars[point],parents,sources[parent_original].triangles[parent])
            except ValueError:pass
            else:valid=[(0.,q)];method='closed-same-original-face-interior-Steiner'
        if not valid:raise ValueError('Unavailable/inconsistent original event ancestry at '+str((stage,point,sorted(parents),method)))
        distance,q=min(valid,key=lambda item:item[0]);result[point]=q
        event=dict(stage=stage,point=q,original_finite_parents=sorted(parents),method=method,replay_m=distance)
        if point in inherited:event['inherited_stage']=inherited[point]['stage']
        events.append(event)
    for point,q in result.items():
        event=next(e for e in events if e['point']==q)
        # No distance welding. If two distinct exact events round to one literal
        # encoded coordinate, representation is ambiguous and the stage stops.
        for key in [point,tuple(float(x)for x in q)]:
            if key in inherited and inherited[key]['point']!=q:raise ValueError('Distinct exact events collide in one Float64 encoding')
            inherited[key]=event
    return result,events
