"""Conform unequal genuine mouth cycles without moving native controls."""
from collections import defaultdict
import numpy as np
from .finite_tissue_boolean import CommonPositionBank
from .shared_triangle_edges import refine_surfaces

def conform_paired_mouth(positions,faces,annulus,loops,normal,tolerance=2e-12):
    v=np.asarray(positions,float);normal=np.asarray(normal,float);normal/=np.linalg.norm(normal);axis=np.eye(3)[np.argmin(abs(normal))];u=axis-normal*(axis@normal);u/=np.linalg.norm(u);w=np.cross(normal,u);basis=np.column_stack([u,w]);center=v[np.concatenate(loops)].mean(0);planarity=float(abs((v[np.concatenate(loops)]-center)@normal).max())
    if planarity>tolerance:raise ValueError('A planar shared construction mouth is required')
    rings=[(v[ring]-center)@basis for ring in loops];angles=[np.mod(np.arctan2(p[:,1],p[:,0]),2*np.pi)for p in rings];common=sorted(np.concatenate(angles));unique=[]
    for angle in common:
        if not unique or angle-unique[-1]>1e-12:unique.append(float(angle))
    if len(unique)>1 and unique[0]+2*np.pi-unique[-1]<=1e-12:unique.pop()
    bank=CommonPositionBank(v,tolerance);plans=defaultdict(set);output=[];maximum=0.
    def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]
    for ids,p,aa in zip(loops,rings,angles):
        ring=[]
        for angle in unique:
            delta=abs(np.angle(np.exp(1j*(aa-angle))));i=int(delta.argmin())
            if delta[i]<=1e-12:ring.append(int(ids[i]));continue
            d=np.array([np.cos(angle),np.sin(angle)]);e=np.roll(p,-1,axis=0)-p;denom=cross(d,e);t=np.divide(cross(p,d),denom,out=np.full(len(p),np.nan),where=abs(denom)>1e-18);radius=np.divide(cross(p,e),denom,out=np.full(len(p),np.nan),where=abs(denom)>1e-18);hit=np.flatnonzero((t>0)&(t<1)&(radius>0))
            if len(hit)!=1:raise ValueError('Mouth must be star-shaped about its common centre')
            i=int(hit[0]);a,b=int(ids[i]),int(ids[(i+1)%len(ids)]);point=v[a]+t[i]*(v[b]-v[a]);control=bank.insert(point);residual=float(np.linalg.norm(point-v[a]-t[i]*(v[b]-v[a])));maximum=max(maximum,residual);plans[tuple(sorted((a,b)))].add(control);ring.append(control)
        if not set(map(int,ids))<=set(ring):raise ValueError('A native mouth control was lost')
        output.append(np.array(ring,int))
    surfaces,parents,count=refine_surfaces(bank,{'wall':np.asarray(faces,int),'annulus':np.asarray(annulus,int)},plans)
    return bank.positions,surfaces['wall'],surfaces['annulus'],output,center,dict(nativeControlsLiteral=np.array_equal(bank.positions[:len(v)],v),sourceCycleControls=[len(q)for q in loops],sharedCycleControls=len(unique),addedEdgeControls=len(bank.positions)-len(v),maximumPlanarityResidualMeters=planarity,maximumCollinearResidualMeters=maximum,refinedSourceFaces=count,wallFaceParents=parents['wall'])
