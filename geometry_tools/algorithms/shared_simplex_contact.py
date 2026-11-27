"""Exact stored-coordinate check that adjacent triangles meet only as declared.

Closed triangle predicates often omit topological neighbours. This separate
gate rejects an overlap hidden by a newly shared edge or vertex.
"""
from fractions import Fraction as F
import numpy as np
from .exact_triangle_separation import exact_segment_triangle_intersection


def _point(p):return tuple(F(float(x))for x in p)
def _sub(a,b):return tuple(x-y for x,y in zip(a,b))
def _cross(a,b):return (a[1]*b[2]-a[2]*b[1],a[2]*b[0]-a[0]*b[2],a[0]*b[1]-a[1]*b[0])

def _coplanar_segment(a,b,triangle):
    a,b=_point(a),_point(b);tri=list(map(_point,triangle))
    normal=_cross(_sub(tri[1],tri[0]),_sub(tri[2],tri[0]));drop=max(range(3),key=lambda k:abs(normal[k]));axes=[k for k in range(3)if k!=drop]
    orient=lambda p,q,r:(q[axes[0]]-p[axes[0]])*(r[axes[1]]-p[axes[1]])-(q[axes[1]]-p[axes[1]])*(r[axes[0]]-p[axes[0]])
    sign=1 if orient(*tri)>0 else-1;low,high=F(0),F(1)
    for x,y in zip(tri,tri[1:]+tri[:1]):
        c0=sign*orient(x,y,a);c1=sign*orient(x,y,b);d=c1-c0
        if d==0:
            if c0<0:return []
        elif d>0:low=max(low,-c0/d)
        else:high=min(high,-c0/d)
        if low>high:return []
    return [tuple(x+t*(y-x)for x,y in zip(a,b))for t in [low,high]]

def only_shared_simplex(first,second):
    a=np.asarray(first,float);b=np.asarray(second,float)
    common=sorted(set(map(_point,a))&set(map(_point,b)))
    if len(common)not in [1,2]:raise ValueError('One explicit shared point or edge is required')
    def allowed(p):
        if len(common)==1:return p==common[0]
        x,y=common;d=_sub(y,x);q=_sub(p,x)
        return not any(_cross(d,q))and F(0)<=sum(u*v for u,v in zip(q,d))<=sum(v*v for v in d)
    for tri,other in [(a,b),(b,a)]:
        for x,y in zip(tri,np.roll(tri,-1,axis=0)):
            hit=exact_segment_triangle_intersection(x,y,other)
            if not hit['intersects']:continue
            if hit['coplanar']:points=_coplanar_segment(x,y,other)
            else:
                t=F(int(hit['parameter_numerator']),int(hit['parameter_denominator']))
                points=[tuple(u+t*(v-u)for u,v in zip(_point(x),_point(y)))]
            if any(not allowed(p)for p in points):
                return dict(passed=False,sharedDimension=len(common)-1,exactStoredCoordinateDecision=True)
    return dict(passed=True,sharedDimension=len(common)-1,exactStoredCoordinateDecision=True)
