"""A bounded minimum-area strip between two directed, nonplanar boundary loops.

The contour-to-grid path construction follows the triangular-tile principle of
Fuchs, Kedem and Uselton (1977), doi:10.1145/359842.359846. This implementation
optimizes one supplied seam and an explicit arclength band; it does not claim
their full toroidal-cycle optimization. Optional forbidden transitions support
independent geometric rejection. A shortest path is NOT an embeddedness proof:
callers must audit the complete finite cell, obstacles and occupied volume.

All supplied points and loop edges are literal. No disk or center is added.
"""
from __future__ import annotations
import numpy as np


class AnnularStripError(ValueError):
    pass


def minimum_area_annular_strip(positions,first,second,*,arclength_band=.25,
                               forbidden_transitions=(),maximum_states=250000):
    v=np.asarray(positions,float);a=np.asarray(first,int);b=np.asarray(second,int)
    if len(a)<3 or len(b)<3 or len(np.unique(a))!=len(a)or len(np.unique(b))!=len(b):
        raise AnnularStripError('Two literal simple directed loops are required')
    if np.intersect1d(a,b).size or not np.isfinite(v).all():
        raise AnnularStripError('Finite disjoint loops are required')
    if not 0<arclength_band<=1:raise AnnularStripError('An explicit positive arclength band ≤1 is required')
    n,m=len(a),len(b)
    if (n+1)*(m+1)>maximum_states:raise AnnularStripError('Finite state budget exceeded')
    pa,pb=v[a],v[b];da=np.roll(pa,-1,axis=0)-pa;db=np.roll(pb,-1,axis=0)-pb
    la,lb=np.linalg.norm(da,axis=1),np.linalg.norm(db,axis=1)
    if np.any(la==0)or np.any(lb==0):raise AnnularStripError('A boundary edge collapsed')
    ta=np.r_[0.,np.cumsum(la)]/la.sum();tb=np.r_[0.,np.cumsum(lb)]/lb.sum()
    ca=np.linalg.norm(np.cross(da[:,None],np.r_[pb,pb[:1]][None]-pa[:,None]),axis=2)/2
    cb=np.linalg.norm(np.cross(np.r_[pa,pa[:1]][:,None]-pb[None],db[None]),axis=2)/2
    ca[ca==0]=np.inf;cb[cb==0]=np.inf
    for direction,i,j in forbidden_transitions:
        if direction==0:ca[i,j]=np.inf
        elif direction==1:cb[i,j]=np.inf
        else:raise AnnularStripError('Unknown grid transition')
    allowed=abs(ta[:,None]-tb[None])<=arclength_band
    costs=np.full((n+1,m+1),np.inf);previous=np.full((n+1,m+1),-1,np.int8);costs[0,0]=0.
    for i in range(n+1):
        for j in range(m+1):
            if not allowed[i,j]or(i==0 and j==0):continue
            ac=costs[i-1,j]+ca[i-1,j]if i else np.inf
            bc=costs[i,j-1]+cb[i,j-1]if j else np.inf
            if ac<=bc and np.isfinite(ac):costs[i,j]=ac;previous[i,j]=0
            elif np.isfinite(bc):costs[i,j]=bc;previous[i,j]=1
    if not np.isfinite(costs[n,m]):raise AnnularStripError('No path in the explicitly bounded seam/band search')
    faces=[];transitions=[];i,j=n,m
    while i or j:
        direction=int(previous[i,j])
        if direction==0:
            faces.append([a[i%n],a[(i-1)%n],b[j%m]]);transitions.append((0,i-1,j));i-=1
        elif direction==1:
            faces.append([a[i%n],b[(j-1)%m],b[j%m]]);transitions.append((1,i,j-1));j-=1
        else:raise AnnularStripError('Broken predecessor chain')
    faces=np.asarray(faces[::-1],int);transitions=np.asarray(transitions[::-1],int)
    assert len(faces)==n+m
    return faces,transitions,dict(method='Fixed-seam minimum-area directed grid strip',
        source='https://dl.acm.org/doi/10.1145/359842.359846',nonplanar_extension_requires_full_geometric_checks=True,
        states=(n+1)*(m+1),reachable_band_states=int(allowed.sum()),arclength_band=float(arclength_band),
        forbidden_transitions=len(forbidden_transitions),triangles=len(faces),area_m2=float(costs[n,m]),
        positions_changed=0,boundary_edges_changed=0,no_lumen_disk=True,anatomical_acceptance=False)
