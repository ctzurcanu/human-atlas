"""Retrangulate coplanar Boolean slivers without moving any control.

Only an interior diagonal may change. The four controls must lie in one
plane within the explicitly supplied construction-error bound. Every new
triangle must have the same orientation in the actual float32 embedding.
"""
import numpy as np

def repair_diagonals(positions,faces,origin,plane_tolerance=2e-12,*,stored_positions=None,corner_normals=None):
    v=np.asarray(positions,float);q=(v-origin).astype('<f4').astype(float)+origin if stored_positions is None else np.asarray(stored_positions,float);f=np.array(faces,int,copy=True);edge={};changes=[]
    cn=None if corner_normals is None else np.array(corner_normals,float,copy=True)
    if cn is not None and(cn.shape!=(len(f),3,3)or not np.isfinite(cn).all()):raise ValueError('Finite face-corner normals required')
    if q.shape!=v.shape or not np.isfinite(q).all():raise ValueError('Explicit finite stored controls must match the source bank')
    def edges(t):return [tuple(sorted((int(a),int(b))))for a,b in zip(t,np.roll(t,-1))]
    def add(i):
        for e in edges(f[i]):edge.setdefault(e,set()).add(i)
    def remove(i):
        for e in edges(f[i]):edge[e].remove(i)
    for i in range(len(f)):add(i)
    for iteration in range(len(f)):
        a=v[f];b=q[f];n=np.cross(a[:,1]-a[:,0],a[:,2]-a[:,0]);qn=np.cross(b[:,1]-b[:,0],b[:,2]-b[:,0]);invalid=np.einsum('ij,ij->i',n,qn)<=0
        if cn is not None:invalid|=np.any(np.einsum('fki,fi->fk',cn,qn)<=0,axis=1)
        bad=np.flatnonzero(invalid)
        if not len(bad):return f,dict(diagonalsChanged=len(changes),maximumPlaneResidualMeters=max([c['planeResidualMeters']for c in changes],default=0.),originalControlsLiteral=True,boundaryEdgesLiteral=True,changes=changes)
        candidates=[]
        for ti in bad:
            for k in range(3):
                x,y,c=map(int,np.roll(f[ti],-k));others=edge[tuple(sorted((x,y)))]-{int(ti)}
                if len(others)!=1:continue
                j=next(iter(others));d=next(int(i)for i in f[j]if i not in[x,y]);newedge=tuple(sorted((c,d)))
                if edge.get(newedge):continue
                big=v[f[j]];plane=np.cross(big[1]-big[0],big[2]-big[0]);length=np.linalg.norm(plane)
                if length==0:continue
                plane/=length;residual=float(abs((v[[x,y,c,d]]-big[0])@plane).max())
                if residual>plane_tolerance:continue
                new=np.array([[c,d,y],[d,c,x]],int);p=q[new];normal=np.cross(p[:,1]-p[:,0],p[:,2]-p[:,0]);projected=normal@plane
                if np.any(projected<=0):continue
                old=f[[ti,j]];sources=[]
                for k,t in enumerate(new):
                    row=[]
                    for vertex in t:
                        source=next((tuple(map(int,s))for s in np.argwhere(old==vertex)if s[0]==k),None)
                        if source is None:source=tuple(map(int,np.argwhere(old==vertex)[0]))
                        row.append(source)
                    sources.append(row)
                if cn is not None:
                    values=np.array([[cn[[ti,j][s[0]],s[1]]for s in row]for row in sources])
                    if np.any(np.einsum('fki,fi->fk',values,normal)<=0):continue
                quality=float(np.min(projected/np.linalg.norm(np.roll(p,-1,axis=1)-p,axis=2).max(1)**2))
                candidates.append((quality,int(ti),j,new,residual,sources))
        if not candidates:raise ValueError(f'{len(bad)} float32 inverted faces have no permitted coplanar diagonal repair')
        _,i,j,new,residual,sources=max(candidates,key=lambda c:c[0]);old=f[[i,j]].copy();remove(i);remove(j);f[i],f[j]=new;add(i);add(j)
        if cn is not None:
            previous=cn[[i,j]].copy();cn[[i,j]]=np.array([[previous[s[0],s[1]]for s in row]for row in sources])
        changes.append(dict(faceIds=[i,j],oldFaces=old.tolist(),newFaces=new.tolist(),planeResidualMeters=residual,newCornerSources=sources))
    raise ValueError('Diagonal repair failed to converge')
