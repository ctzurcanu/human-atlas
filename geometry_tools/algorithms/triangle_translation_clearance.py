"""Finite exact-support candidates for translating a triangle past another.

This is a geometric separating-plane certificate, not a route-side choice
or anatomical interface classifier. All source corners enter each bound.
"""
import numpy as np

def separating_translation(first,second,direction,*,clearance_m=0.):
    a,b,d=map(lambda v:np.asarray(v,float),(first,second,direction))
    if a.shape!=(3,3) or b.shape!=(3,3) or d.shape!=(3,) or not all(np.isfinite(x).all() for x in [a,b,d]):raise ValueError('Two finite triangles and one direction required')
    if not np.isfinite(clearance_m) or clearance_m<0:raise ValueError('Finite nonnegative clearance required')
    length=float(np.linalg.norm(d))
    if length==0:raise ValueError('Nonzero translation direction required')
    world_scale=float(max(np.abs(a).max(),np.abs(b).max()))
    d=d/length;origin=(a.mean(0)+b.mean(0))*.5;a=a-origin;b=b-origin
    ea=np.roll(a,-1,axis=0)-a;eb=np.roll(b,-1,axis=0)-b
    na=np.cross(ea[0],ea[1]);nb=np.cross(eb[0],eb[1])
    if np.linalg.norm(na)==0 or np.linalg.norm(nb)==0:raise ValueError('Degenerate source triangles need explicit treatment')
    # Triangle normals and edge-cross-edge axes cover ordinary intersections;
    # in-plane edge normals also cover coplanar faces. Direction itself is an
    # additional always-valid support candidate, not a necessary SAT axis.
    axes=np.vstack([na,nb,np.cross(ea[:,None],eb[None]).reshape(-1,3),np.cross(na,ea),np.cross(nb,eb),d])
    norms=np.linalg.norm(axes,axis=1);axes=axes[norms>0]/norms[norms>0,None]
    slopes=axes@d;axes[slopes<0]*=-1;slopes=np.abs(slopes)
    valid=slopes>128*np.finfo(float).eps;axes=axes[valid];slopes=slopes[valid]
    lower=(a@axes.T).min(0);upper=(b@axes.T).max(0)
    required=np.maximum((upper-lower+clearance_m)/slopes,0)
    i=int(np.argmin(required));axis=axes[i];motion=float(required[i])
    gap=float(((a+motion*d)@axis).min()-(b@axis).max())
    rounding=512*np.finfo(float).eps*(world_scale+np.abs(origin).max()+np.abs(a).max()+np.abs(b).max()+motion+clearance_m)
    if gap+rounding<clearance_m:raise ValueError('Source support certificate lost to arithmetic')
    return motion,dict(complete=True,direction=d.tolist(),separatingAxis=axis.tolist(),axisSlope=float(slopes[i]),minimumCommonTranslation_m=motion,certifiedProjectedGap_m=gap,requestedClearance_m=float(clearance_m),supportArithmeticBound_m=float(rounding),candidateAxes=len(axes),allTriangleCornersIncluded=True,validForLargerNonuniformCornerTranslationsAlongDirection=True,anatomicalAcceptance=False)
