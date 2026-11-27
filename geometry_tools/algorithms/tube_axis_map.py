"""Coherent section-frame transport near a tube axis, with a literal fixed tail."""
import numpy as np
from scipy.spatial.transform import Rotation,Slerp
from .tube_sections import arclength,rotation_minimizing_frames,unit

def axis_transport(points,old_axis,new_axis,rotations):
    p=np.asarray(points,float);old=np.asarray(old_axis,float);new=np.asarray(new_axis,float);r=np.asarray(rotations,float)
    if old.shape!=new.shape or old.shape[1:]!=(3,)or r.shape!=(len(old),3,3):raise ValueError('Corresponding axes and proper section rotations required')
    if not np.allclose(np.linalg.det(r),1.,rtol=0,atol=1e-12)or not np.allclose(np.swapaxes(r,1,2)@r,np.eye(3),rtol=0,atol=1e-12):raise ValueError('Reflection, shear or singular section frame')
    a=old[:-1];d=np.diff(old,axis=0);length2=np.sum(d*d,axis=1);arclength(old)
    t=np.clip(np.einsum('pki,ki->pk',p[:,None]-a,d)/length2,0,1);near=a[None]+t[:,:,None]*d[None]
    k=np.argmin(np.sum((p[:,None]-near)**2,axis=2),axis=1);q=t[np.arange(len(p)),k]
    start=near[np.arange(len(p)),k];end=new[k]+q[:,None]*(new[k+1]-new[k])
    rotations_at=Slerp(np.arange(len(old)),Rotation.from_matrix(r))(k+q).as_matrix()
    out=end+np.einsum('pij,pj->pi',rotations_at,p-start)
    fixed=np.all(old[k]==new[k],axis=1)&np.all(old[k+1]==new[k+1],axis=1)&np.all(r[k]==np.eye(3),axis=(1,2))&np.all(r[k+1]==np.eye(3),axis=(1,2))
    out[fixed]=p[fixed]
    return out,rotations_at

def endpoint_course(old_axis,prefix_points,prefix_tangents,last_station):
    old=np.asarray(old_axis,float);prefix=np.asarray(prefix_points,float);t=unit(prefix_tangents);k=int(last_station)
    if k<1 or k>=len(old)or not np.array_equal(prefix[-1],old[k]):raise ValueError('Prefix must terminate at the actual fixed source station')
    from .tube_sections import path_tangents
    old_t=path_tangents(old);before=rotation_minimizing_frames(old,old_t)
    arc=arclength(prefix);sample=arclength(old[:k+1]);sample=sample/sample[-1]*arc[-1]
    new=np.column_stack([np.interp(sample,arc,prefix[:,j])for j in range(3)]);new[[0,-1]]=prefix[[0,-1]]
    new_t=unit(np.column_stack([np.interp(sample,arc,t[:,j])for j in range(3)]));new_t[-1]=old_t[k]
    after=rotation_minimizing_frames(new,new_t,initial_reference=before[0,:,0],final_reference=before[k,:,0])
    rotations=np.broadcast_to(np.eye(3),(len(old),3,3)).copy();rotations[:k+1]=after@np.swapaxes(before[:k+1],1,2)
    rotations[k]=np.eye(3);new=np.vstack([new,old[k+1:]])
    return new,rotations
