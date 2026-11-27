"""Proper-pose plane intervals and analytic local derivatives.

Float64 replay, not exact arithmetic. A separating plane is sufficient only
for its recorded pair; complete actual triangle contacts remain required.
"""
import numpy as np
from scipy import sparse
from scipy.spatial.transform import Rotation
from .rigid_section_twists import replay_section_twists
from .original_singular_pose_rows import skew,right_rotation_jacobian

class ExactPosePairPlanes:
    def __init__(self,first_points,second_points,first_groups,second_groups,axes,minimum_interval_m,centers,levers):
        self.a,self.b,self.axis=[np.asarray(v,float) for v in [first_points,second_points,axes]]
        self.ga,self.gb=[np.asarray(v) for v in [first_groups,second_groups]]
        self.minimum=np.asarray(minimum_interval_m,float);self.centers=np.asarray(centers,float);self.levers=np.asarray(levers,float);self.count=len(self.centers)
        if self.a.ndim!=2 or self.a.shape[1]!=3 or self.b.shape!=self.a.shape or self.axis.shape!=self.a.shape or self.minimum.shape!=(len(self.a),) or self.ga.shape!=self.minimum.shape or self.gb.shape!=self.minimum.shape or self.centers.shape!=(len(self.levers),3):
            raise ValueError('Explicit paired points, plane ownership and proper-pose cells required')
        if any(g.dtype.kind not in 'iu' or np.any(g< -1) or np.any(g>=self.count) for g in [self.ga,self.gb]) or np.any(self.levers<=0) or np.any(self.minimum<0) or not all(np.isfinite(v).all() for v in [self.a,self.b,self.axis,self.minimum,self.centers,self.levers]) or np.any(abs(np.linalg.norm(self.axis,axis=1)-1)>1e-10):
            raise ValueError('Finite source data, unit planes and fixed-cell sentinel required')
    def values(self,parameters):
        a,_=replay_section_twists(self.a,self.ga,self.centers,self.levers,parameters)
        b,_=replay_section_twists(self.b,self.gb,self.centers,self.levers,parameters)
        return np.einsum('ij,ij->i',self.axis,a-b)
    def violations(self,parameters,tolerance_m=1e-10):
        values=self.values(parameters);bad=np.flatnonzero(self.minimum-values>tolerance_m)
        return bad,dict(complete=True,checkedExactConstraints=len(values),violatingExactConstraints=len(bad),maximumIntervalDeficit_m=float(np.maximum(self.minimum-values,0).max(initial=0)),actualProperRotationsReplayed=True,originalPlanesAndTargetsUnchanged=True)
    def local_rows(self,parameters,ids):
        x=np.asarray(parameters,float);ids=np.asarray(ids,int);phi=x[self.count:]/self.levers[:,None]
        rotations=Rotation.from_rotvec(phi).as_matrix();right=right_rotation_jacobian(phi);rows=[];cols=[];data=[]
        for sign,p,groups in [(1.,self.a,self.ga),(-1.,self.b,self.gb)]:
            group=groups[ids];live=group>=0;g=group[live];row=np.repeat(np.flatnonzero(live),3)
            v=p[ids[live]]-self.centers[g];J=-rotations[g]@skew(v)@right[g]/self.levers[g,None,None]
            translation=sign*self.axis[ids[live]];rotation=np.einsum('ni,nij->nj',translation,J)
            rows.extend([row,row]);cols.extend([(g[:,None]*3+np.arange(3)).ravel(),(3*self.count+g[:,None]*3+np.arange(3)).ravel()]);data.extend([translation.ravel(),rotation.ravel()])
        G=sparse.coo_matrix((np.concatenate(data),(np.concatenate(rows),np.concatenate(cols))),shape=(len(ids),6*self.count)).tocsr();G.eliminate_zeros()
        norm=np.sqrt(np.asarray(G.power(2).sum(1)).ravel());values=self.values(x)[ids]
        if np.any((norm<1e-12)&(self.minimum[ids]-values>1e-10)):raise ValueError('A violated interval has no resolved pose derivative')
        scale=np.maximum(1.,1/np.maximum(norm,1e-12));rhs=(self.minimum[ids]-values+G@x.ravel())*scale
        return sparse.diags(scale)@G,rhs,dict(complete=True,localLinearizationOnly=True,exactFinitePoseReplayRequired=True)
