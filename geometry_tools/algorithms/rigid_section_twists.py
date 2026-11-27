"""Linearized section poses with exact proper-rotation replay.

One six-parameter pose acts on every control of each declared rigid cell.
Shared cells may span explicitly coupled source mouths; proximity creates no
membership. A linear solve is only a proposal: replay uses Rodrigues rotations
and its actual surfaces must pass the original geometric acceptance gates.
"""
import numpy as np
from scipy import sparse
from scipy.spatial.transform import Rotation


def section_twist_lift(points, groups, *, group_count, parameter_radius_m,
                       maximum_angle_radians):
    p=np.asarray(points,float);g=np.asarray(groups,int)
    if p.ndim!=2 or p.shape[1]!=3 or g.shape!=(len(p),) or not np.isfinite(p).all():
        raise ValueError('Finite source controls and explicit cell membership required')
    if group_count<1 or np.any(g<0) or np.any(g>=group_count) or set(g)!=set(range(group_count)):
        raise ValueError('Every declared cell must have source controls')
    if not 0<parameter_radius_m<=.01 or not 0<maximum_angle_radians<=np.pi/6:
        raise ValueError('Explicit finite pose bounds required')
    centers=np.array([np.unique(p[g==i],axis=0).mean(0) for i in range(group_count)])
    radii=np.array([np.linalg.norm(p[g==i]-centers[i],axis=1).max() for i in range(group_count)])
    # Each translation and scaled rotation has its own parameter ball. This
    # lever bounds exact rotation motion to <= 20% of that ball radius, and
    # the angle to the independently declared angular cap.
    levers=np.maximum(5*radii,parameter_radius_m/maximum_angle_radians)
    r=(p-centers[g])/levers[g,None]
    skew=np.zeros((len(p),3,3))
    skew[:,0,1]=-r[:,2];skew[:,0,2]=r[:,1]
    skew[:,1,0]=r[:,2];skew[:,1,2]=-r[:,0]
    skew[:,2,0]=-r[:,1];skew[:,2,1]=r[:,0]
    row=np.arange(3*len(p)).reshape(-1,3)
    rr=np.r_[row.ravel(),np.broadcast_to(row[:,:,None],skew.shape).ravel()]
    cc=np.r_[(g[:,None]*3+np.arange(3)).ravel(),
             (3*group_count+g[:,None,None]*3+np.arange(3)[None,None,:]+np.zeros((len(p),3,1),int)).ravel()]
    values=np.r_[np.ones(3*len(p)),-skew.ravel()]
    lift=sparse.coo_matrix((values,(rr,cc)),shape=(3*len(p),6*group_count)).tocsr()
    lift.eliminate_zeros()
    return lift,centers,levers,dict(complete=True,groupCount=group_count,
        maximumExactControlMotionBound_m=1.2*parameter_radius_m,
        parameterRadius_m=parameter_radius_m,maximumAngleRadians=maximum_angle_radians,
        proximityMembershipInferred=False,exactReplayRequired=True)


def replay_section_twists(points, groups, centers, levers, parameters):
    p=np.asarray(points,float);g=np.asarray(groups,int);centers=np.asarray(centers,float)
    n=len(centers);x=np.asarray(parameters,float)
    if x.shape!=(2*n,3) or g.shape!=(len(p),) or np.any(g<-1) or np.any(g>=n):
        raise ValueError('Corresponding bounded section parameters required')
    rotation=Rotation.from_rotvec(x[n:]/np.asarray(levers)[:,None]).as_matrix()
    live=g>=0;out=p.copy();s=g[live]
    out[live]=np.einsum('nij,nj->ni',rotation[s],p[live]-centers[s])+centers[s]+x[s]
    # Fixed raw controls are never round-tripped through any floating pose.
    assert np.array_equal(out[~live],p[~live])
    return out,dict(complete=True,properRotations=True,fixedControlsLiteral=True,
        maximumAngleRadians=float(np.linalg.norm(x[n:]/np.asarray(levers)[:,None],axis=1).max()),
        maximumActualMotion_m=float(np.linalg.norm(out-p,axis=1).max()),
        exactProfileDistancesPreserved=True,anatomicalAcceptance=False)
