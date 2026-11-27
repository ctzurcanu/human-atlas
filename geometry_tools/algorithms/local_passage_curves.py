"""Compare predeclared local transverse loops; infer no anatomical interface."""
import numpy as np
from scipy.optimize import linear_sum_assignment


def _segments(curves, component):
    ids=set(component['controls']);p=np.array(curves['points_m'])
    edges=np.array([e for e in curves['edges']if e[0]in ids and e[1]in ids],int)
    return p[edges]


def _directed_bound(first,second,spacing):
    samples=[];maximum_step=0.
    for a,b in first:
        length=float(np.linalg.norm(b-a));count=max(1,int(np.ceil(length/spacing)))
        samples.append(a+np.linspace(0,1,count+1)[:,None]*(b-a));maximum_step=max(maximum_step,length/count)
    p=np.concatenate(samples);a=second[:,0];v=second[:,1]-a;denom=np.einsum('ij,ij->i',v,v)
    if np.any(denom==0):raise ValueError('A transverse curve cannot contain zero-length edges')
    t=np.clip(np.einsum('psi,si->ps',p[:,None]-a,v)/denom,0,1)
    distance=np.linalg.norm(p[:,None]-(a+t[:,:,None]*v),axis=2).min(1)
    # Distance to a closed set is 1-Lipschitz; midpoint interval radius bounds
    # the unsampled intervals, with a conservative Float64 arithmetic margin.
    return float(distance.max()+maximum_step/2+1e-9)


def compare_local_passage_curves(before,after,*,maximum_displacement_m,
                                 maximum_length_ratio=1.25,spacing_m=.00005):
    old=before['components'];new=after['components'];failures=[]
    if any(not c['closedUnbranched']for c in old+new):failures.append('Open or branched intersection curve')
    if len(new)>len(old):failures.append('Additional intersection loops')
    if failures:return dict(passed=False,failures=failures,matched=[],removedOldLoops=[])
    if not new:return dict(passed=True,failures=[],matched=[],removedOldLoops=list(range(len(old))))
    costs=np.full((len(new),len(old)),1e6);evidence={}
    for i,n in enumerate(new):
        ns=_segments(after,n)
        for j,o in enumerate(old):
            os=_segments(before,o);ratio=n['length_m']/o['length_m']
            if ratio>maximum_length_ratio:continue
            if np.linalg.norm(np.array(n['centroid_m'])-o['centroid_m'])>2*maximum_displacement_m:continue
            bound=max(_directed_bound(ns,os,spacing_m),_directed_bound(os,ns,spacing_m))
            if bound<=maximum_displacement_m:costs[i,j]=bound;evidence[i,j]=dict(newLoop=i,oldLoop=j,hausdorffUpperBound_m=bound,lengthRatio=ratio)
    first,second=linear_sum_assignment(costs);matched=[]
    for i,j in zip(first,second):
        if (int(i),int(j))not in evidence:failures.append('New loop cannot match a predeclared local passage')
        else:matched.append(evidence[int(i),int(j)])
    return dict(passed=not failures,failures=failures,matched=matched,
                removedOldLoops=sorted(set(range(len(old)))-set(map(int,second))),
                anatomicalInterfaceEstablished=False)
