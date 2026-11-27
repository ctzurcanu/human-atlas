"""Bounded local rigid interpolation of explicit corresponding source controls."""
import time
import numpy as np
from scipy.spatial import cKDTree


def interpolate_rigid_field(points, source, target, *, neighbours=48, support_m=None,
                            wall_seconds=15., batch_size=512):
    start=time.monotonic();x,p,q=[np.asarray(v,float)for v in (points,source,target)]
    if any(v.ndim!=2 or v.shape[1]!=3 or not np.isfinite(v).all()for v in (x,p,q))or p.shape!=q.shape or len(p)<4:
        raise ValueError('Finite query points and at least four corresponding 3D source controls required')
    if not 4<=neighbours<=128 or not 0<wall_seconds<=20 or not 1<=batch_size<=2048:
        raise ValueError('Explicit finite neighbourhood, memory batch and native deadline required')
    unique,rep,inv=np.unique(p,axis=0,return_index=True,return_inverse=True)
    if not np.array_equal(q,q[rep][inv]):
        raise ValueError('One literal source control has contradictory target poses')
    p,q=unique,q[rep]
    if support_m is not None and not (len(support_m)==2 and 0<support_m[0]<support_m[1]):
        raise ValueError('An explicit ordered inner/outer support radius is required')
    tree=cKDTree(p);result=x.copy();exact_count=0;max_motion=0.
    for first in range(0,len(x),batch_size):
        if time.monotonic()-start>wall_seconds:raise TimeoutError('Local rigid field native deadline')
        part=x[first:first+batch_size];distance,ids=tree.query(part,k=min(neighbours,len(p)))
        exact=distance[:,0]==0;exact_count+=int(exact.sum());posed=part.copy();posed[exact]=q[ids[exact,0]]
        active=~exact
        if support_m is not None:active &= distance[:,0]<support_m[1]
        if np.any(active):
            d=distance[active];idx=ids[active]
            weights=(d[:,0,None]/d)**4;weights/=weights.sum(1)[:,None]
            old,new=p[idx],q[idx];a=np.einsum('bn,bni->bi',weights,old);b=np.einsum('bn,bni->bi',weights,new)
            H=np.einsum('bni,bnj,bn->bij',old-a[:,None],new-b[:,None],weights)
            u,s,vt=np.linalg.svd(H);v=vt.transpose(0,2,1);ut=u.transpose(0,2,1)
            signs=np.ones((len(v),3));signs[:,-1]=np.linalg.det(v@ut)
            R=(v*signs[:,None])@ut
            if np.any(np.linalg.det(R)<.999999):raise ValueError('Local rigid interpolation reflected source geometry')
            transformed=np.einsum('bij,bj->bi',R,part[active]-a)+b
            if support_m is not None:
                t=np.clip((support_m[1]-d[:,0])/(support_m[1]-support_m[0]),0,1)
                fade=t*t*(3-2*t);transformed=part[active]+fade[:,None]*(transformed-part[active])
            posed[active]=transformed
        result[first:first+len(part)]=posed
        max_motion=max(max_motion,float(np.linalg.norm(posed-part,axis=1).max(initial=0)))
    return result,dict(complete=True,sourceControls=len(p),queryPoints=len(x),exactInterpolatedControls=exact_count,
        neighbours=neighbours,batchSize=batch_size,support_m=support_m,maximumMotion_m=max_motion,
        localRotationsProper=True,globalInjectivityProven=False,anatomicalAcceptance=False,
        elapsedSeconds=time.monotonic()-start)
