"""Proper section transport with exact prescribed branch-end orientations."""
import numpy as np
from .tube_sections import path_tangents, rotation_minimizing_frames


def branch_anchored_frames(source_centres, target_centres, start_rotation, end_rotation):
    """Transport source section obliquity and spread endpoint roll by arclength.

    The source rings need not be perpendicular to their course. Returned proper
    matrices act on their original offsets. This does not prove tube validity.
    """
    p,q=np.asarray(source_centres,float),np.asarray(target_centres,float)
    ends=np.asarray([start_rotation,end_rotation],float)
    if p.shape!=q.shape or p.ndim!=2 or p.shape[1]!=3 or len(p)<2:
        raise ValueError('Corresponding ordered source and target chains required')
    if ends.shape!=(2,3,3) or not np.isfinite(ends).all() or not np.allclose(ends.transpose(0,2,1)@ends,np.eye(3),atol=1e-11,rtol=0) or np.any(np.linalg.det(ends)<1-1e-11):
        raise ValueError('Proper branch attachment orientations required')
    if np.array_equal(p,q) and np.array_equal(ends,np.repeat(np.eye(3)[None],2,axis=0)):
        return np.repeat(np.eye(3)[None],len(p),axis=0)
    old_t=path_tangents(p);new_t=path_tangents(q)
    new_t[0]=ends[0]@old_t[0];new_t[-1]=ends[1]@old_t[-1]
    old=rotation_minimizing_frames(p,old_t)
    new=rotation_minimizing_frames(q,new_t,initial_reference=ends[0]@old[0,:,0],
                                  final_reference=ends[1]@old[-1,:,0])
    R=new@old.transpose(0,2,1)
    if not np.allclose(R[[0,-1]],ends,atol=1e-11,rtol=0):
        raise ValueError('Transport could not retain the branch orientations')
    R[0],R[-1]=ends
    return R
