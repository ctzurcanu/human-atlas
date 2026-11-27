"""Bounded diagonal changes on near-planar, unprotected Boolean slivers.

Controls and the outer patch boundary stay literal. Each changed patch has an
explicit source-planarity bound. Full encoded topology/crossing gates follow.
"""
import numpy as np

def regularize_transport_slivers(source,stored,faces,fixed,*,quality_floor=1e-5,plane_tolerance_m=1e-7,max_changes=512,allowed_faces=None):
    v=np.asarray(source,float);q=np.asarray(stored,float);f=np.array(faces,int,copy=True);edge={};changes=[]
    def keys(t):return [tuple(sorted((int(a),int(b))))for a,b in zip(t,np.roll(t,-1))]
    def add(i):
        for e in keys(f[i]):edge.setdefault(e,set()).add(i)
    def remove(i):
        for e in keys(f[i]):edge[e].remove(i)
    def quality(t):
        return np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/np.linalg.norm(t-np.roll(t,-1,axis=1),axis=2).max(1)**2
    allowed=set(map(int,allowed_faces))if allowed_faces is not None else None
    for i in range(len(f)):add(i)
    for iteration in range(max_changes):
        score=quality(q[f]);bad=np.flatnonzero(score<quality_floor);candidates=[]
        if allowed is not None:bad=np.array([i for i in bad if i in allowed],int)
        for ti in bad:
            for k in range(3):
                x,y,c=map(int,np.roll(f[ti],-k));others=edge[tuple(sorted((x,y)))]-{int(ti)}
                if len(others)!=1:continue
                j=next(iter(others));d=next(int(i)for i in f[j]if i not in[x,y]);quad=[x,y,c,d]
                if np.any(np.asarray(fixed)[quad])or edge.get(tuple(sorted((c,d)))):continue
                big=v[f[j]];plane=np.cross(big[1]-big[0],big[2]-big[0]);length=np.linalg.norm(plane)
                if length==0:continue
                plane/=length;residual=float(abs((v[quad]-big[0])@plane).max())
                if residual>plane_tolerance_m:continue
                new=np.array([[c,d,y],[d,c,x]],int);p=q[new];normal=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);big=q[f[j]];reference=np.cross(big[1]-big[0],big[2]-big[0])
                if np.any(normal@reference<=0):continue
                improved=float(quality(p).min())
                if improved<=min(score[ti],score[j])*2:continue
                candidates.append((improved,int(ti),j,new,residual))
        if not candidates:return f,dict(changes=changes,sourceControlsLiteral=True,patchBoundaryLiteral=True,maximumSourcePlanarityMetres=max((c['sourcePlanarityMetres']for c in changes),default=0.))
        _,i,j,new,residual=max(candidates,key=lambda c:c[0]);old=f[[i,j]].copy();remove(i);remove(j);f[i],f[j]=new;add(i);add(j)
        if allowed is not None:allowed.update([i,j])
        changes.append(dict(faceIds=[i,j],oldFaces=old.tolist(),newFaces=new.tolist(),sourcePlanarityMetres=residual))
    raise ValueError('Sliver diagonal-change budget reached')
