"""Recover occupied components from oriented closed Boolean boundary shells.

A Boolean engine may return an outer positive shell and separate negative
cavity shells from ``decompose``. Those shells are not separate occupied
volumes. This helper assigns each negative shell to its smallest enclosing
positive shell and classifies query points against that signed region.
Input shells must already be nonintersecting, closed oriented boundaries.
"""
import numpy as np
import igl


def occupied_shell_regions(shells,queries):
    q=np.asarray(queries,float)
    if q.ndim!=2 or q.shape[1]!=3 or not np.isfinite(q).all():raise ValueError('Finite [N,3] queries required')
    volumes=[]
    for mesh in shells:
        tri=mesh.triangles;origin=mesh.vertices.mean(0);t=tri-origin
        volume=float(np.einsum('ij,ij->i',t[:,0],np.cross(t[:,1],t[:,2])).sum()/6)
        if volume==0:raise ValueError('Zero-volume shell cannot define a solid component')
        volumes.append(volume)
    positive=[i for i,v in enumerate(volumes)if v>0];negative=[i for i,v in enumerate(volumes)if v<0];children={i:[]for i in positive}
    def winding(index,points):
        mesh=shells[index];w=igl.winding_number(mesh.vertices,mesh.faces,points)
        if np.any(np.minimum(abs(w),abs(abs(w)-1))>1e-7):raise ValueError('A query or shell witness is on/ambiguous against another boundary')
        return w
    for i in negative:
        witness=shells[i].triangles[0,0][None]
        owners=[j for j in positive if abs(winding(j,witness)[0])>.5]
        if not owners:raise ValueError('A negative boundary has no enclosing positive shell')
        children[min(owners,key=lambda j:volumes[j])].append(i)
    rows=[]
    for i in positive:
        inside=winding(i,q)>.5
        for child in children[i]:inside&=abs(winding(child,q))<.5
        volume=volumes[i]+sum(volumes[j]for j in children[i])
        if volume<=0:raise ValueError('Cavity boundaries exhaust their positive component')
        rows.append(dict(outer_shell=i,negative_shells=children[i],volume_m3=volume,query_indices=np.flatnonzero(inside).tolist()))
    return rows
