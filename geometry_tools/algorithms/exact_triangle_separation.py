"""Exact separating-axis decision on literal nondegenerate Float64 triangles.

Each binary float is converted to an integer on the pair's common power-of-two
grid. Separation signs have no distance epsilon or cancellation threshold.
The result concerns stored surface geometry, never anatomical acceptance.
"""
import math
import numpy as np


def exact_triangle_separation(first,second):
    aa=np.asarray(first,float);bb=np.asarray(second,float)
    if aa.shape!=(3,3)or bb.shape!=(3,3)or not np.isfinite(np.r_[aa.ravel(),bb.ravel()]).all():
        raise ValueError('Two finite literal triangles are required')
    ratios=[float(x).as_integer_ratio()for x in np.r_[aa.ravel(),bb.ravel()]]
    denominator=max(d for _,d in ratios)
    values=[n*(denominator//d)for n,d in ratios]
    a=[tuple(values[i:i+3])for i in range(0,9,3)];b=[tuple(values[i:i+3])for i in range(9,18,3)]
    sub=lambda x,y:tuple(u-v for u,v in zip(x,y))
    cross=lambda x,y:(x[1]*y[2]-x[2]*y[1],x[2]*y[0]-x[0]*y[2],x[0]*y[1]-x[1]*y[0])
    dot=lambda x,y:sum(u*v for u,v in zip(x,y))
    ea=[sub(a[(i+1)%3],a[i])for i in range(3)];eb=[sub(b[(i+1)%3],b[i])for i in range(3)]
    na=cross(ea[0],ea[1]);nb=cross(eb[0],eb[1])
    if not any(na)or not any(nb):raise ValueError('Degenerate triangles require their own dimensional predicate')
    axes=[na,nb]+[cross(x,y)for x in ea for y in eb]+[cross(na,e)for e in ea]+[cross(nb,e)for e in eb]
    for index,axis in enumerate(axes):
        if not any(axis):continue
        pa=[dot(p,axis)for p in a];pb=[dot(p,axis)for p in b]
        gap=max(min(pb)-max(pa),min(pa)-max(pb))
        if gap>0:
            distance=float(gap)/math.sqrt(float(dot(axis,axis)))/denominator
            return dict(disjoint=True,separating_axis_index=index,integer_gap=str(gap),binary_grid_denominator=str(denominator),
                projected_separation_m=float(distance),exact_stored_coordinate_decision=True)
    return dict(disjoint=False,exact_stored_coordinate_decision=True,interpretation='Closed triangles intersect or touch')


def exact_segment_triangle_intersection(first,last,triangle):
    p=np.asarray([first,last],float);t=np.asarray(triangle,float)
    if p.shape!=(2,3)or t.shape!=(3,3)or not np.isfinite(np.r_[p.ravel(),t.ravel()]).all():
        raise ValueError('A finite segment and literal triangle are required')
    ratios=[float(x).as_integer_ratio()for x in np.r_[p.ravel(),t.ravel()]];denominator=max(d for _,d in ratios)
    values=[n*(denominator//d)for n,d in ratios]
    pts=[tuple(values[i:i+3])for i in range(0,15,3)];a,b,x,y,z=pts
    sub=lambda u,v:tuple(q-r for q,r in zip(u,v))
    cross=lambda u,v:(u[1]*v[2]-u[2]*v[1],u[2]*v[0]-u[0]*v[2],u[0]*v[1]-u[1]*v[0])
    dot=lambda u,v:sum(q*r for q,r in zip(u,v))
    e1=sub(y,x);e2=sub(z,x);normal=cross(e1,e2)
    if not any(normal):raise ValueError('A nondegenerate source triangle is required')
    direction=sub(b,a);h=cross(direction,e2);det=dot(e1,h)
    if det==0:
        if dot(normal,sub(a,x))!=0:return dict(intersects=False,coplanar=False,exact_stored_coordinate_decision=True)
        drop=max(range(3),key=lambda i:abs(normal[i]));keep=[i for i in range(3)if i!=drop]
        project=lambda u:tuple(u[i]for i in keep)
        a,b,x,y,z=map(project,[a,b,x,y,z])
        orient=lambda u,v,w:(v[0]-u[0])*(w[1]-u[1])-(v[1]-u[1])*(w[0]-u[0])
        def inside(q):
            s=[orient(u,v,q)for u,v in [(x,y),(y,z),(z,x)]]
            return all(v>=0 for v in s)or all(v<=0 for v in s)
        hit=inside(a)or inside(b)
        for u,v in [(x,y),(y,z),(z,x)]:
            bbox=all(max(min(a[i],b[i]),min(u[i],v[i]))<=min(max(a[i],b[i]),max(u[i],v[i]))for i in range(2))
            hit=hit or (bbox and orient(a,b,u)*orient(a,b,v)<=0 and orient(u,v,a)*orient(u,v,b)<=0)
        return dict(intersects=bool(hit),coplanar=True,exact_stored_coordinate_decision=True)
    s=sub(a,x);u=dot(s,h);q=cross(s,e1);v=dot(direction,q);parameter=dot(e2,q)
    if det<0:det,u,v,parameter=-det,-u,-v,-parameter
    hit=0<=u<=det and 0<=v<=det and u+v<=det and 0<=parameter<=det
    return dict(intersects=bool(hit),coplanar=False,at_first_endpoint=bool(hit and parameter==0),
        at_last_endpoint=bool(hit and parameter==det),parameter_numerator=str(parameter)if hit else None,
        parameter_denominator=str(det)if hit else None,exact_stored_coordinate_decision=True)
