"""Improve a shared internal surface without moving its controls or rim.

Callers must keep its two tissue owners reciprocal and must audit the final
occupied volumes. Near-planar edge flips are construction candidates, not an
anatomical classifier or a substitute for surface-intersection checks.
"""
from collections import defaultdict
import time
import numpy as np
from .tube_sections import TubeGeometryError


def improve_triangles(source,encoded,faces,*,plane_tolerance_m=1e-6,
                      poor_quality=.01,maximum_passes=8,time_limit=5.,
                      protected_edges=()):
    p=np.asarray(source,float);q=np.asarray(encoded,float);f=np.array(faces,np.int64,copy=True)
    if p.shape!=q.shape or p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all()or not np.isfinite(q).all():
        raise TubeGeometryError('Corresponding finite source and encoded controls required')
    if f.ndim!=2 or f.shape[1]!=3 or not len(f):raise TubeGeometryError('An oriented triangle patch is required')
    if plane_tolerance_m<0 or not np.isfinite(plane_tolerance_m)or maximum_passes<1 or time_limit<=0:
        raise TubeGeometryError('Finite nonnegative deviation and positive iteration/deadline limits required')
    def quality(tri):
        xyz=q[tri];n=np.cross(xyz[...,1,:]-xyz[...,0,:],xyz[...,2,:]-xyz[...,0,:])
        den=np.sum((xyz-np.roll(xyz,-1,axis=-2))**2,axis=(-2,-1))
        return np.divide(2*np.sqrt(3)*np.linalg.norm(n,axis=-1),den,out=np.zeros_like(den),where=den>0)
    def edges(face):return [tuple(sorted((int(a),int(b))))for a,b in zip(face,np.roll(face,-1))]
    owners=defaultdict(set)
    for i,tri in enumerate(f):
        for edge in edges(tri):owners[edge].add(i)
    if any(len(v)>2 for v in owners.values()):raise TubeGeometryError('A patch edge is nonmanifold')
    protected={tuple(sorted(map(int,edge)))for edge in protected_edges}
    if any(len(edge)!=2 or edge not in owners for edge in protected):
        raise TubeGeometryError('Protected edges must belong to the declared patch')
    before_boundary={edge for edge,rows in owners.items()if len(rows)==1}
    start=time.monotonic();flips=[];max_deviation=0.;completed=0
    for iteration in range(maximum_passes):
        if time.monotonic()-start>=time_limit:raise TimeoutError('Shared surface quality improvement reached its deadline')
        quality0=quality(f)
        src=p[f];dst=q[f]
        sn=np.cross(src[:,1]-src[:,0],src[:,2]-src[:,0]);dn=np.cross(dst[:,1]-dst[:,0],dst[:,2]-dst[:,0])
        wrong=np.einsum('ij,ij->i',sn,dn)<=0
        quality0[wrong]=-quality0[wrong]
        candidates=sorted((min(quality0[list(rows)]),edge)for edge,rows in owners.items()
                          if len(rows)==2 and edge not in protected and min(quality0[list(rows)])<poor_quality)
        used=set();changed=0
        for _,edge in candidates:
            if time.monotonic()-start>=time_limit:raise TimeoutError('Shared surface quality improvement reached its deadline')
            rows=owners.get(edge,set())
            if len(rows)!=2 or used.intersection(rows):continue
            i,j=sorted(rows);a,b=edge
            tri=f[i];location=next((k for k in range(3)if set((int(tri[k]),int(tri[(k+1)%3])))==set(edge)),None)
            if location is None:continue
            a,b,c=map(int,(tri[location],tri[(location+1)%3],tri[(location+2)%3]))
            other=f[j];pos=next((k for k in range(3)if int(other[k])==b and int(other[(k+1)%3])==a),None)
            if pos is None:raise TubeGeometryError('A shared patch edge has inconsistent winding')
            d=int(other[(pos+2)%3]);diagonal=tuple(sorted((c,d)))
            if c==d or owners.get(diagonal):continue
            old=f[[i,j]];new=np.array([[c,d,b],[d,c,a]],np.int64)
            original=p[old];normal=np.cross(original[:,1]-original[:,0],original[:,2]-original[:,0]).sum(0)
            size=np.linalg.norm(normal)
            if size==0:continue
            normal/=size;quad=q[[a,b,c,d]];deviation=float(np.abs((quad-quad.mean(0))@normal).max())
            if deviation>plane_tolerance_m:continue
            xyz=q[new];cross=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
            if np.any(cross@normal<=0):continue
            xyz=p[new];cross_source=np.cross(xyz[:,1]-xyz[:,0],xyz[:,2]-xyz[:,0])
            if np.any(cross_source@normal<=0):continue
            old_score=float(quality0[[i,j]].min());new_score=float(quality(new).min())
            if new_score<=old_score*(1+1e-6)+1e-12:continue
            for row in (i,j):
                for key in edges(f[row]):
                    owners[key].remove(row)
                    if not owners[key]:del owners[key]
            f[i],f[j]=new
            for row in (i,j):
                for key in edges(f[row]):owners[key].add(row)
            used.update((i,j));changed+=1;max_deviation=max(max_deviation,deviation)
            flips.append(dict(rows=[i,j],old_diagonal=[a,b],new_diagonal=[c,d],minimum_quality_before=old_score,minimum_quality_after=new_score))
        completed=iteration+1
        if not changed:break
    boundary={edge for edge,rows in owners.items()if len(rows)==1}
    if boundary!=before_boundary:raise TubeGeometryError('The patch boundary changed')
    if any(edge not in owners for edge in protected):raise TubeGeometryError('A protected feature edge changed')
    if any(len(rows)>2 for rows in owners.values()):raise TubeGeometryError('A flip made a patch edge nonmanifold')
    return f,dict(source_coordinates_moved=False,encoded_coordinates_moved=False,
                  patch_boundary_edges_literal=True,vertices_before=len(np.unique(faces)),vertices_after=len(np.unique(f)),
                  protected_feature_edges=len(protected),protected_edges_literal=True,
                  flips=flips,passes=completed,plane_tolerance_m=plane_tolerance_m,
                  maximum_quad_plane_deviation_m=max_deviation,elapsed_seconds=time.monotonic()-start,
                  full_mesh_acceptance=False)
