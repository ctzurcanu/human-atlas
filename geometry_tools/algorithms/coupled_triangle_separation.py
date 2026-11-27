"""Minimum squared corner translations for explicit triangle support axes.

Both triangles may move; protected corners stay literal. This constructs
geometric proposals, not anatomical interfaces or accepted mesh repairs.
"""
import numpy as np


def separating_pair_proposals(vertices, first_ids, second_ids, *, locked=(), gap_m=2e-7):
    p=np.asarray(vertices,float); ia=np.asarray(first_ids,int);ib=np.asarray(second_ids,int)
    if ia.shape!=(3,) or ib.shape!=(3,) or len(set(map(int,np.r_[ia,ib])))!=6:
        raise ValueError('Two disjoint explicit source triangles required')
    if not np.isfinite(p).all() or not np.isfinite(gap_m) or gap_m<=0:
        raise ValueError('Finite source and positive support margin required')
    ids=np.r_[ia,ib];origin=p[ids].mean(0);a=p[ia]-origin;b=p[ib]-origin
    ea=np.roll(a,-1,axis=0)-a;eb=np.roll(b,-1,axis=0)-b
    na=np.cross(ea[0],ea[1]);nb=np.cross(eb[0],eb[1])
    if min(np.linalg.norm(na),np.linalg.norm(nb))<=1e-24:
        raise ValueError('Nondegenerate explicit source facets required')
    axes=np.vstack([na,nb,np.cross(ea[:,None],eb[None]).reshape(-1,3),np.cross(na,ea),np.cross(nb,eb)])
    lengths=np.linalg.norm(axes,axis=1);axes=axes[lengths>0]/lengths[lengths>0,None]
    axes=np.vstack([axes,-axes]);fixed=np.isin(ids,list(locked));rows=[]
    for index,n in enumerate(axes):
        s=np.r_[a@n,b@n];low=float(s.min()-gap_m);high=float(s.max()+gap_m)
        if fixed[:3].any():high=min(high,float(s[:3][fixed[:3]].min()-gap_m*.5))
        if fixed[3:].any():low=max(low,float(s[3:][fixed[3:]].max()+gap_m*.5))
        if low>high:continue
        for _ in range(64):
            middle=(low+high)*.5
            da=np.maximum(middle+gap_m*.5-s[:3],0);db=np.minimum(middle-gap_m*.5-s[3:],0)
            derivative=float(da[~fixed[:3]].sum()+db[~fixed[3:]].sum())
            if derivative>0:high=middle
            else:low=middle
        middle=(low+high)*.5
        delta=np.r_[np.maximum(middle+gap_m*.5-s[:3],0),np.minimum(middle-gap_m*.5-s[3:],0)]
        delta[fixed]=0;values=p[ids]+delta[:,None]*n
        certificate=float(((values[:3]-origin)@n).min()-((values[3:]-origin)@n).max())
        rounding=256*np.finfo(float).eps*max(float(np.abs(p[ids]).max()),gap_m)
        if certificate+rounding<gap_m:continue
        changed=np.flatnonzero(delta!=0)
        rows.append(dict(ids=ids[changed],positions=values[changed],axis=n,
            squaredMotion=float(delta@delta),maximumMotion_m=float(np.abs(delta).max()),
            projectedGap_m=certificate,axisIndex=index,arithmeticBound_m=rounding))
    return sorted(rows,key=lambda r:(r['squaredMotion'],r['maximumMotion_m'],r['axisIndex']))
