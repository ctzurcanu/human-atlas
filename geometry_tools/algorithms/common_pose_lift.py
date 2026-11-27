"""First derivative of explicitly coupled proper poses, preserving source offsets."""
import numpy as np
from scipy import sparse

def common_pose_lift(old_centers,old_levers,new_centers,new_levers,old_to_new):
    old,new=np.asarray(old_centers,float),np.asarray(new_centers,float)
    a,b=np.asarray(old_levers,float),np.asarray(new_levers,float)
    mapping=np.asarray(old_to_new)
    if old.shape!=(len(a),3) or new.shape!=(len(b),3) or mapping.shape!=(len(a),) or mapping.dtype.kind not in 'iu':
        raise ValueError('Explicit matching proper-pose cells required')
    if min(len(a),len(b))<1 or mapping.min()<0 or mapping.max()>=len(b) or np.any(a<=0) or np.any(b<=0) or not all(np.isfinite(x).all() for x in [old,new,a,b]):
        raise ValueError('Finite positive levers and complete source-cell ownership required')
    rows=[];cols=[];values=[];n=len(a);m=len(b)
    def block(r,c,M):
        for i,j in zip(*np.nonzero(M)):
            rows.append(r+i);cols.append(c+j);values.append(M[i,j])
    for g,h in enumerate(mapping):
        delta=old[g]-new[h];x,y,z=delta
        cross=np.array([[0.,-z,y],[z,0.,-x],[-y,x,0.]])
        block(3*g,3*h,np.eye(3))
        block(3*g,3*(m+h),-cross/b[h])
        block(3*(n+g),3*(m+h),np.eye(3)*a[g]/b[h])
    return sparse.coo_matrix((values,(rows,cols)),shape=(6*n,6*m)).tocsr()
