"""Local plane subdivision with shared-edge propagation and corner ancestry.

Only requested facets create plane intersections. Adjacent unrequested facets
receive any required edge subdivision, keeping the entire tissue complex
conforming without extending the local cage throughout the body.
"""
import numpy as np


def refine_selected_plane(positions,faces,*,normal,coordinate,origin=(0.,0.,0.),
                          active_faces=None,plane_tolerance_m=1e-12):
    p=np.asarray(positions,float);f=np.asarray(faces,int);n=np.asarray(normal,float)
    o=np.asarray(origin,float);coordinate=float(coordinate);tol=float(plane_tolerance_m)
    if p.ndim!=2 or p.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or n.shape!=(3,)or o.shape!=(3,):
        raise ValueError('Physical 3D controls, triangles and an explicit plane are required')
    if not np.isfinite(np.r_[p.ravel(),n,o,coordinate,tol]).all()or np.linalg.norm(n)==0 or tol<0:
        raise ValueError('Finite plane data and nonnegative classification tolerance required')
    if len(f)and(f.min()<0 or f.max()>=len(p)):raise ValueError('Invalid triangle indices')
    active=np.ones(len(f),bool)if active_faces is None else np.asarray(active_faces,bool)
    if active.shape!=(len(f),):raise ValueError('One active flag per source facet required')
    length=np.linalg.norm(n);n=n/length;coordinate/=length
    d=np.asarray((p.astype(np.longdouble)-o.astype(np.longdouble))@n.astype(np.longdouble)-np.longdouble(coordinate),float)
    d[np.abs(d)<=tol]=0
    face_d=d[f]
    edges=np.stack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]],axis=1)
    crossing=(d[edges[:,:,0]]*d[edges[:,:,1]]<0)&active[:,None]
    selected=np.unique(np.sort(edges[crossing],axis=1),axis=0)
    identity=np.eye(3)
    if not len(selected):
        return dict(positions=p.copy(),faces=f.copy(),parent_faces=np.arange(len(f)),
            corner_barycentrics=np.broadcast_to(identity,(len(f),3,3)).copy(),
            evidence=dict(new_controls=0,added_faces=0,selected_faces_split=0,propagated_faces=0,
                          original_controls_literal=True))
    a,b=selected.T;parameter=d[a]/(d[a]-d[b]);new=p[a]+parameter[:,None]*(p[b]-p[a])
    vertices=np.r_[p,new];base=len(p)
    if base>2_000_000_000:raise ValueError('Explicit edge-key size limit exceeded')
    keys=selected[:,0]*base+selected[:,1];ordered=np.sort(edges,axis=2)
    edge_keys=ordered[:,:,0]*base+ordered[:,:,1]
    at=np.searchsorted(keys,edge_keys);take=(at<len(keys))&(keys[np.minimum(at,len(keys)-1)]==edge_keys)
    cuts=np.where(take,base+at,-1);count=take.sum(1)
    if np.any(count>2):raise ValueError('A triangle cannot have three strictly crossing plane edges')
    zero=np.flatnonzero(count==0);one=np.flatnonzero(count==1);two=np.flatnonzero(count==2)
    k=np.argmax(take[one],axis=1);order=(k[:,None]+np.arange(3))%3
    fo=np.take_along_axis(f[one],order,axis=1);da=np.take_along_axis(face_d[one],order,axis=1)
    aa,bb,cc=fo.T;ab=cuts[one,k]
    first=np.stack([np.c_[aa,ab,cc],np.c_[ab,bb,cc]],axis=1).reshape(-1,3)
    ba,bc,bd=identity[order[:,0]],identity[order[:,1]],identity[order[:,2]]
    t=da[:,0]/(da[:,0]-da[:,1]);mid=ba+t[:,None]*(bc-ba)
    first_bary=np.stack([np.stack([ba,mid,bd],axis=1),np.stack([mid,bc,bd],axis=1)],axis=1).reshape(-1,3,3)
    dt=face_d[two];lone=np.where(np.sum(dt>0,axis=1)==1,np.argmax(dt,axis=1),np.argmin(dt,axis=1))
    order=(lone[:,None]+np.arange(3))%3;ft=np.take_along_axis(f[two],order,axis=1)
    dd=np.take_along_axis(dt,order,axis=1);aa,bb,cc=ft.T
    ab=cuts[two,lone];ac=cuts[two,(lone+2)%3]
    if np.any(ab<0)or np.any(ac<0):raise ValueError('Two cut edges must meet at the plane outlier')
    second=np.stack([np.c_[aa,ab,ac],np.c_[ab,bb,cc],np.c_[ab,cc,ac]],axis=1).reshape(-1,3)
    ba,bc,bd=identity[order[:,0]],identity[order[:,1]],identity[order[:,2]]
    tab=dd[:,0]/(dd[:,0]-dd[:,1]);tac=dd[:,0]/(dd[:,0]-dd[:,2])
    pab=ba+tab[:,None]*(bc-ba);pac=ba+tac[:,None]*(bd-ba)
    second_bary=np.stack([np.stack([ba,pab,pac],axis=1),np.stack([pab,bc,bd],axis=1),
                         np.stack([pab,bd,pac],axis=1)],axis=1).reshape(-1,3,3)
    result=np.r_[f[zero],first,second];parents=np.r_[zero,np.repeat(one,2),np.repeat(two,3)]
    bary=np.concatenate([np.broadcast_to(identity,(len(zero),3,3)),first_bary,second_bary])
    source=p[f[parents]];target=vertices[result]
    before=np.cross(source[:,1]-source[:,0],source[:,2]-source[:,0])
    after=np.cross(target[:,1]-target[:,0],target[:,2]-target[:,0])
    if np.any(np.einsum('ni,ni->n',before,after)<=0):
        raise ValueError('Selective plane subdivision produced a collapsed or reversed facet')
    changed=count>0
    return dict(positions=vertices,faces=result,parent_faces=parents,corner_barycentrics=bary,
        evidence=dict(new_controls=len(new),added_faces=len(result)-len(f),
            selected_faces_split=int(np.sum(changed&active)),propagated_faces=int(np.sum(changed&~active)),
            original_controls_literal=True,shared_edge_intersections_stored_once=True,
            source_surface_unchanged=True,plane_classification_tolerance_m=tol,
            storage_and_anatomical_review_pending=True))
