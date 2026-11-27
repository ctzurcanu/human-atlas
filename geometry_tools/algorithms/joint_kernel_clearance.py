"""Small shared translation fields fitted to simultaneous clearance rows.

Signed gaps and their increasing-distance normals come from explicit measured
surfaces. The optimizer produces a bounded proposal, never accepted anatomy.
Actual surface queries, coarse triangles and source interfaces must be replayed.
"""
import time
import numpy as np
from scipy.optimize import minimize
from scipy.linalg import null_space
from scipy.spatial import cKDTree
from scipy.spatial.distance import cdist


class GaussianTranslationBasis:
    def __init__(self,centres,width_m,*,frame,free_below,fixed_above,background=1e-4):
        self.centres=np.array(centres,float,copy=True)
        if self.centres.ndim!=2 or self.centres.shape[1]!=3 or not 4<=len(self.centres)<=192 or not np.isfinite(self.centres).all():
            raise ValueError('Between four and 192 finite explicit field centres required')
        if not np.isfinite([width_m,free_below,fixed_above,background]).all()or not .005<=width_m<=.03 or not free_below<fixed_above or not 1e-8<=background<=.01 or not isinstance(frame,str)or not frame.strip():
            raise ValueError('Explicit bounded field width, root support, background and frame required')
        self.width=float(width_m);self.frame=frame;self.free_below=float(free_below);self.fixed_above=float(fixed_above);self.background=float(background)

    def weights(self,points,*,gradients=False):
        p=np.asarray(points,float)
        if p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all():raise ValueError('Finite query points required')
        K=np.exp(-cdist(p,self.centres,'sqeuclidean')/(2*self.width**2));W=K/(K.sum(1)[:,None]+self.background)
        u=np.clip((self.fixed_above-p[:,1])/(self.fixed_above-self.free_below),0,1);fade=u*u*(3-2*u)
        if not gradients:return W*fade[:,None]
        delta=self.centres[None]-p[:,None];mean=np.einsum('nk,nki->ni',W,delta)
        G=W[:,:,None]*(delta-mean[:,None])/self.width**2
        df=-6*u*(1-u)/(self.fixed_above-self.free_below);df[(u<=0)|(u>=1)]=0
        G*=fade[:,None,None];G[:,:,1]+=W*df[:,None]
        return W*fade[:,None],G

    def map_and_jacobians(self,points,coefficients):
        c=np.asarray(coefficients,float)
        if c.shape!=self.centres.shape or not np.isfinite(c).all():raise ValueError('One finite translation per field centre required')
        W,G=self.weights(points,gradients=True)
        return np.asarray(points)+W@c,np.eye(3)[None]+np.einsum('nki,kj->nji',G,c)

    def evidence(self,coefficients):
        return dict(type='normalized-gaussian-translation',centres_m=self.centres.tolist(),coefficients_m=np.asarray(coefficients).tolist(),width_m=self.width,background=self.background,freeBelow_m=self.free_below,fixedAbove_m=self.fixed_above,frame=self.frame,continuousInjectivityCertified=False,anatomicalAcceptance=False)


def propose_joint_step(points,basis,gaps,normals,row_point_ids,*,requested_gap_m=.00005,
                       maximum_step_m=.0006,maximum_iterations=80,wall_seconds=4.,stationary_point_ids=()):
    """Fit both surfaces together with a strong no-deepening penalty.

    Unknown surface sides must be omitted by the caller. A true nonlinear
    surface replay still decides whether the proposal can be retained.
    Translation coefficients are normalized after the bounded convex solve so
    every query displacement is at most ``maximum_step_m``.
    """
    start=time.monotonic();p=np.asarray(points,float);gap=np.asarray(gaps,float);n=np.asarray(normals,float);ids=np.asarray(row_point_ids,int)
    if gap.ndim!=1 or n.shape!=(len(gap),3)or ids.shape!=gap.shape or np.any(ids<0)or np.any(ids>=len(p)):
        raise ValueError('Explicit gap/normal/point correspondence required')
    if not np.isfinite(gap).all()or not np.isfinite(n).all()or not np.allclose(np.linalg.norm(n,axis=1),1,atol=1e-8):
        raise ValueError('Only known finite gaps and unit increasing-clearance normals allowed')
    if not 0<requested_gap_m<=.001 or not 0<maximum_step_m<=.002 or not 1<=maximum_iterations<=160 or not 0<wall_seconds<=10:
        raise ValueError('Finite explicit gap, step and solver limits required')
    W=basis.weights(p);B=W[ids];soft=(requested_gap_m-gap)/maximum_step_m
    hard=(np.minimum(gap,0)-gap)/maximum_step_m
    k=len(basis.centres);near=cKDTree(basis.centres).query(basis.centres,k=min(5,k))[1][:,1:]
    edges=np.unique(np.sort(np.c_[np.repeat(np.arange(k),near.shape[1]),near.ravel()],axis=1),axis=0)
    protected=np.unique(np.asarray(stationary_point_ids,int))
    if np.any(protected<0)or np.any(protected>=len(p)):raise ValueError('Stationary controls must be actual query indices')
    Z=np.eye(k);constraints=np.empty((0,k))
    if len(protected):
        pw,pg=basis.weights(p[protected],gradients=True)
        constraints=np.r_[pw,pg.transpose(0,2,1).reshape(-1,k)]
        length=np.linalg.norm(constraints,axis=1);constraints=constraints[length>1e-14]/length[length>1e-14,None]
        if len(constraints):Z=null_space(constraints,rcond=1e-11)
    reduced=B@Z
    if not Z.shape[1]:return np.zeros((k,3)),dict(complete=True,proposalOnly=True,optimizerConverged=True,optimizerMessage='No nonstationary field degrees of freedom',stationaryControls=len(protected),elapsedSeconds=time.monotonic()-start,anatomicalAcceptance=False)
    calls=0
    def objective(x):
        nonlocal calls
        if time.monotonic()-start>wall_seconds:raise TimeoutError('Joint clearance solver deadline; no proposal returned')
        calls+=1;z=x.reshape(Z.shape[1],3);c=Z@z;v=reduced@z;projection=np.einsum('ni,ni->n',v,n)
        a=np.maximum(soft-projection,0);h=np.maximum(hard-projection,0)
        residual=a+1000*h
        value=.5*float(a@a+1000*(h@h))/len(gap)+.0001*float((c*c).sum())/k
        gradient=-(reduced.T@(residual[:,None]*n))/len(gap)
        difference=c[edges[:,0]]-c[edges[:,1]];value+=.001*float((difference*difference).sum())/len(edges)
        g=.002*difference/len(edges);regularization=.0002*c/k;np.add.at(regularization,edges[:,0],g);np.add.at(regularization,edges[:,1],-g)
        return value,(gradient+Z.T@regularization).ravel()
    size=Z.shape[1]*3
    result=minimize(objective,np.zeros(size),jac=True,method='L-BFGS-B',bounds=[(-1.,1.)]*size,options=dict(maxiter=maximum_iterations,maxfun=maximum_iterations*3,ftol=1e-12,gtol=1e-8,maxls=20))
    coefficients=Z@result.x.reshape(Z.shape[1],3);maximum=float(np.linalg.norm(coefficients,axis=1).max());coefficients*=maximum_step_m/max(maximum,1.)
    q,J=basis.map_and_jacobians(p,coefficients)
    assert np.linalg.norm(q-p,axis=1).max()<=maximum_step_m*(1+1e-9)
    if len(protected):
        np.testing.assert_allclose(q[protected],p[protected],atol=1e-13,rtol=0)
        np.testing.assert_allclose(J[protected],np.broadcast_to(np.eye(3),J[protected].shape),atol=1e-10,rtol=0)
    return coefficients,dict(complete=True,proposalOnly=True,optimizerConverged=bool(result.success),optimizerMessage=str(result.message),iterations=int(result.nit),objectiveCalls=calls,stationaryControls=len(protected),fieldDegreesOfFreedom=int(size),maximumRequestedStep_m=maximum_step_m,maximumActualStep_m=float(np.linalg.norm(q-p,axis=1).max()),elapsedSeconds=time.monotonic()-start,anatomicalAcceptance=False)
