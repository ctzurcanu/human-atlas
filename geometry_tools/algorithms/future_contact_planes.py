"""Source-bound separating planes from explicitly rejected triangle trials.

The original pose supplies a separating plane; the rejected trial identifies
which of its inequalities must enter the next linear solve. Negative right
hand sides retain available clearance. Shared source cells receive zero gap,
so this helper never opens an attached source junction to invent separation.
"""
import numpy as np


def future_contact_rows(vertices,faces,owner,feedback,protected_cells,*,gap_m,complete_pair_planes=False):
    v=np.asarray(vertices,float);f=np.asarray(faces,int);o=np.asarray(owner,int)
    for key,value in [('vertices',v),('faces',f),('owner',o)]:
        if not np.array_equal(feedback[key],value):raise ValueError('Future contact feedback belongs to another '+key)
    pairs=np.asarray(feedback['pairs'],int);shift=np.asarray(feedback['probe_shift'],float)
    if pairs.ndim!=2 or pairs.shape[1]!=2 or shift.shape!=(int(o.max())+1,3)or not np.isfinite(shift).all():
        raise ValueError('Actual source face pairs and a finite cell trial required')
    if len(pairs)and(pairs.min()<0 or pairs.max()>=len(f)):raise ValueError('Feedback face outside source array')
    if not len(pairs):return np.empty((0,2),int),np.empty((0,3)),np.empty(0),dict(pairs=0,rows=0)
    a,b=v[f[pairs[:,0]]],v[f[pairs[:,1]]];oa,ob=o[f[pairs[:,0]]],o[f[pairs[:,1]]]
    ea=a[:,[1,2,0]]-a;eb=b[:,[1,2,0]]-b
    axes=np.concatenate([np.cross(ea[:,0],ea[:,1])[:,None],np.cross(eb[:,0],eb[:,1])[:,None],
        np.cross(ea[:,:,None],eb[:,None]).reshape(-1,9,3)],axis=1)
    size=np.linalg.norm(axes,axis=2);valid=size>1e-20
    axes=np.divide(axes,size[:,:,None],out=np.zeros_like(axes),where=valid[:,:,None])
    axes=np.concatenate([axes,-axes],axis=1);valid=np.concatenate([valid,valid],axis=1)
    origin=a[:,0]
    def gaps(x,y):
        px=np.einsum('nki,nji->nkj',x-origin[:,None],axes);py=np.einsum('nki,nji->nkj',y-origin[:,None],axes)
        return py.min(1)-px.max(1)
    original=gaps(a,b);probe=gaps(a+shift[oa],b+shift[ob]);valid&=original>=-1e-10
    supported=valid.any(1);score=np.where(valid,probe,-np.inf);chosen=score.argmax(1)
    normal=axes[np.arange(len(pairs)),chosen]
    shared=np.any(oa[:,:,None]==ob[:,None,:],axis=(1,2))
    available=np.maximum(original[np.arange(len(pairs)),chosen],0.)
    gap=np.where(shared,0.,np.minimum(gap_m,.5*available))
    groups=[];directions=[];values=[];fixed=np.asarray(protected_cells,bool)
    for i in range(3):
        for j in range(3):
            first,second=oa[:,i],ob[:,j]
            rhs=np.minimum(np.einsum('ni,ni->n',a[:,i]-b[:,j],normal)+gap,0.)
            predicted=np.einsum('ni,ni->n',shift[second]-shift[first],normal)
            # A subsequent solve can move a previously safe corner. Retaining
            # every corner inequality makes this a complete separating-plane
            # constraint, rather than just repairing the old probe's active
            # corners. The sparse legacy mode remains explicit for callers
            # whose saved contracts use those individual contact cuts.
            active=np.ones(len(pairs),bool)if complete_pair_planes else rhs>predicted+1e-9
            keep=supported&active&(first!=second)&(~fixed[first]|~fixed[second])
            groups.append(np.c_[first[keep],second[keep]]);directions.append(normal[keep]);values.append(rhs[keep])
    rows=np.concatenate(groups);normals=np.concatenate(directions);rhs=np.concatenate(values)
    return rows,normals,rhs,dict(pairs=len(pairs),rows=len(rows),completePairPlanes=bool(complete_pair_planes),sourcePoseSeparatingPlanes=int(supported.sum()),
        unsupportedPairs=int((~supported).sum()),sharedSourceCellPairs=int(shared.sum()),
        sharedSourceCellGap_m=0.,negativeClearanceBoundsRetained=int(np.sum(rhs<0)),anatomicalAcceptance=False)
