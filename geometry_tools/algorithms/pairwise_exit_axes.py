"""Full triangle-interval exits along separately declared local route directions."""
import numpy as np

def pairwise_exit_axes(first,second,direction,*,clearance_m=0.):
    """SAT escape planes with least common translation along an explicit ray.

    For each literal triangle pair, chooses a full-interval separating plane
    that can be reached in the positive supplied direction. Subsequent variable
    section transport must enforce all nine vertex inequalities, not assume the
    triangles move rigidly together. This is not a whole-mesh collision proof.
    """
    a,b=np.asarray(first,float),np.asarray(second,float);direction=np.asarray(direction,float)
    if a.ndim!=3 or a.shape[1:]!=(3,3)or b.shape!=a.shape or direction.shape!=(len(a),3)or not np.isfinite(a).all()or not np.isfinite(b).all()or not np.isfinite(direction).all()or np.any(np.linalg.norm(direction,axis=1)<1e-12) or not np.isfinite(clearance_m)or clearance_m<0:
        raise ValueError('Finite matching triangles and an explicit positive escape ray required')
    direction=direction/np.linalg.norm(direction,axis=1)[:,None]
    ea=np.roll(a,-1,axis=1)-a;eb=np.roll(b,-1,axis=1)-b;na=np.cross(ea[:,0],ea[:,1]);nb=np.cross(eb[:,0],eb[:,1])
    axes=np.concatenate([na[:,None],nb[:,None],np.cross(ea[:,:,None],eb[:,None]).reshape(-1,9,3),np.cross(na[:,None],ea),np.cross(nb[:,None],eb)],axis=1)
    lengths=np.linalg.norm(axes,axis=2);axes=np.divide(axes,lengths[:,:,None],out=np.zeros_like(axes),where=lengths[:,:,None]>1e-25)
    rate=np.einsum('pai,pi->pa',axes,direction);axes*=np.where(rate<0,-1.,1.)[:,:,None];rate=abs(rate)
    pa=np.einsum('pvc,pac->pav',a,axes);pb=np.einsum('pvc,pac->pav',b,axes);gap=pa.min(2)-pb.max(2)
    amount=np.divide(np.maximum(clearance_m-gap,0.),rate,out=np.full_like(rate,np.inf),where=rate>1e-10)
    chosen=amount.argmin(1);index=np.arange(len(a))
    if not np.isfinite(amount[index,chosen]).all():raise ValueError('No finite escape plane along the declared ray')
    return axes[index,chosen],amount[index,chosen]


