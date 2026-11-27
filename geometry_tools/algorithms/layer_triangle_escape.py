"""Directed triangle exit planes within an explicitly supplied tissue layer.

Projected triangle-edge normals provide the 2D Minkowski separating axes.
The least displacement along the supplied in-layer ray selects a plane; every
vertex inequality still has to be enforced and actual poses independently
checked. This chooses geometric planes, not an anatomical vessel route.
"""
import numpy as np


def layer_triangle_exit_axes(first,second,layer_normals,escape_directions,*,clearance_m=0.):
    a,b,n,d=[np.asarray(x,float).copy() for x in [first,second,layer_normals,escape_directions]]
    if a.ndim!=3 or a.shape[1:]!=(3,3) or b.shape!=a.shape or n.shape!=(len(a),3) or d.shape!=n.shape:
        raise ValueError('Explicit paired triangles, layer normals and directed in-layer rays required')
    if not all(np.isfinite(x).all() for x in [a,b,n,d]) or not np.isfinite(clearance_m) or clearance_m<0:
        raise ValueError('Finite source geometry and nonnegative gap required')
    lengths=np.linalg.norm(n,axis=1)
    if np.any(lengths<1e-12):raise ValueError('Resolved layer normals required')
    n/=lengths[:,None]
    d-=np.einsum('ij,ij->i',d,n)[:,None]*n
    lengths=np.linalg.norm(d,axis=1)
    if np.any(lengths<1e-8):raise ValueError('Escape direction must resolve a direction in the layer')
    d/=lengths[:,None]
    edges=np.concatenate([np.roll(a,-1,axis=1)-a,np.roll(b,-1,axis=1)-b],axis=1)
    axes=np.concatenate([np.cross(n[:,None],edges),d[:,None]],axis=1)
    lengths=np.linalg.norm(axes,axis=2)
    axes=np.divide(axes,lengths[:,:,None],out=np.zeros_like(axes),where=lengths[:,:,None]>1e-20)
    rate=np.einsum('pki,pi->pk',axes,d)
    axes*=np.where(rate<0,-1.,1.)[:,:,None];rate=abs(rate)
    pa=np.einsum('pvi,pki->pkv',a,axes);pb=np.einsum('pvi,pki->pkv',b,axes)
    gap=pa.min(2)-pb.max(2)
    amount=np.divide(np.maximum(clearance_m-gap,0),rate,out=np.full_like(rate,np.inf),where=rate>1e-10)
    chosen=amount.argmin(1);row=np.arange(len(a))
    if not np.isfinite(amount[row,chosen]).all():raise ValueError('No finite directed in-layer exit plane')
    result=axes[row,chosen]
    return result,amount[row,chosen],dict(complete=True,pairs=len(a),
        maximumLayerNormalComponent=float(abs(np.einsum('ij,ij->i',result,n)).max(initial=0)),
        method='Projected edge-normal Minkowski planes, minimum directed in-layer interval exit',
        sourceGeometryModified=False,actualVertexAndPoseVerificationRequired=True,anatomicalRouteEstablished=False)
