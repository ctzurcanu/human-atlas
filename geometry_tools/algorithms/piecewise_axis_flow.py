"""Globally invertible piecewise affine axis maps on conformed triangles.

The primary coordinate is strictly monotone. The other coordinates receive
explicit piecewise linear shears. A triangle must first be split at every
crossed knot; applying the field to unrefined triangles is rejected.
"""
import numpy as np


def conformed_axis_flow(positions,faces,*,axis,knots,targets):
    p=np.asarray(positions,float);f=np.asarray(faces,int);k=int(axis)
    x=np.asarray(knots,float);t=np.asarray(targets,float)
    if p.ndim!=2 or p.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or k not in (0,1,2) \
            or x.ndim!=1 or len(x)<2 or t.shape!=(len(x),3):
        raise ValueError('Explicit axis knots and three target components are required')
    if not np.isfinite(np.r_[p.ravel(),x,t.ravel()]).all() or np.any(np.diff(x)<=0):
        raise ValueError('Finite, strictly increasing source knots are required')
    slopes=np.diff(t,axis=0)/np.diff(x)[:,None]
    if np.any(slopes[:,k]<=0):raise ValueError('Primary coordinate must be strictly increasing')
    if len(f) and (f.min()<0 or f.max()>=len(p)):raise ValueError('Invalid triangle indices')
    q=p[f,k]
    crossed=np.searchsorted(x,q.max(1),side='left')>np.searchsorted(x,q.min(1),side='right')
    if np.any(crossed):raise ValueError('Source triangles must conform to every crossed cage knot')
    region=np.clip(np.searchsorted(x,p[:,k],side='right')-1,0,len(x)-2)
    value=t[region]+(p[:,k]-x[region])[:,None]*slopes[region]
    out=p+value;out[:,k]=value[:,k]
    face_region=np.clip(np.searchsorted(x,q.mean(1),side='right')-1,0,len(x)-2)
    return dict(positions=out,face_slopes=slopes[face_region],
        evidence=dict(continuous_global_homeomorphism=True,
            analytic_determinant_lower_bound=float(slopes[:,k].min()),
            source_triangles_conform_to_all_knots=True,face_count_unchanged=True,
            float_storage_and_anatomical_acceptance_pending=True))


def transport_axis_corner_normals(normals,face_slopes,*,axis):
    n=np.asarray(normals,float);s=np.asarray(face_slopes,float);k=int(axis)
    if n.ndim!=3 or n.shape[1:]!=(3,3) or s.shape!=(len(n),3) or k not in (0,1,2):
        raise ValueError('Corresponding corner normals and positive affine face slopes required')
    if not np.isfinite(np.r_[n.ravel(),s.ravel()]).all() or np.any(s[:,k]<=0):
        raise ValueError('Finite, positive affine face slopes required')
    out=n.copy();other=[j for j in range(3)if j!=k]
    out[:,:,k]=(n[:,:,k]-np.einsum('ti,tki->tk',s[:,other],n[:,:,other]))/s[:,None,k]
    return out
