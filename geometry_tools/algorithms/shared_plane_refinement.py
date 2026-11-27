"""Conform shared physical triangle charts to one explicit cage plane.

This changes tessellation only. Every crossed physical edge receives exactly
one cached intersection; output corner barycentrics retain field ancestry.
It adds no tissue, cap or anatomical attachment. Call once per saved stage.
"""
import numpy as np


def safe_plane_coordinate(values,desired,*,search_radius_m=.00075,minimum_separation_m=1e-10):
    """Choose an explicit nearby plane without creating near-vertex slivers.

    A desired coordinate in a broad empty interval or outside the mesh is
    already safe, even when no interval midpoint lies in the search window.
    """
    q=np.unique(np.asarray(values,float));desired=float(desired)
    radius=float(search_radius_m);guard=float(minimum_separation_m)
    if not len(q) or not np.isfinite(np.r_[q,desired,radius,guard]).all() or radius<0 or guard<0:
        raise ValueError('Finite coordinates and nonnegative plane-selection limits required')
    gaps=np.diff(q);mid=(q[:-1]+q[1:])*.5
    take=np.flatnonzero((np.abs(mid-desired)<=radius)&(gaps>=2*guard))
    if len(take):return float(mid[take[np.argmax(gaps[take])]])
    if np.min(np.abs(q-desired))>=guard:return desired
    raise ValueError('No plane satisfies the requested search window and separation guard')


def lift_refined_positions(source_positions,source_faces,refinement):
    """Lift cage-space subdivisions through literal physical corner ancestry.

    This avoids round-tripping original physical controls through a rotated
    cage frame. Shared new controls use one cached parent-corner witness.
    """
    p=np.asarray(source_positions,float);f=np.asarray(source_faces,int)
    faces=np.asarray(refinement['faces'],int);parents=np.asarray(refinement['parent_faces'],int)
    bary=np.asarray(refinement['corner_barycentrics'],float)
    out=np.empty_like(np.asarray(refinement['positions'],float));out[:len(p)]=p
    flat=faces.ravel();new=flat>=len(p)
    ids,first=np.unique(flat[new],return_index=True)
    if not np.array_equal(ids,np.arange(len(p),len(out))):
        raise ValueError('Every appended cage control requires an actual source-face witness')
    corner=np.flatnonzero(new)[first]
    out[ids]=np.einsum('ni,nij->nj',bary.reshape(-1,3)[corner],p[f[parents[corner//3]]])
    return out


def refine_at_axis_plane(positions,faces,*,axis,coordinate,minimum_vertex_distance_m=1e-10):
    p=np.asarray(positions,float);f=np.asarray(faces,int);axis=int(axis);plane=float(coordinate)
    if p.ndim!=2 or p.shape[1]!=3 or f.ndim!=2 or f.shape[1]!=3 or axis not in [0,1,2]:
        raise ValueError('Three-dimensional positions, triangles and a valid axis are required')
    if not np.isfinite(np.r_[p.ravel(),plane,minimum_vertex_distance_m]).all() or minimum_vertex_distance_m<0:
        raise ValueError('Finite coordinates and a nonnegative plane separation are required')
    if len(f) and (f.min()<0 or f.max()>=len(p)):raise ValueError('Invalid source indices')
    distance=p[:,axis]-plane
    if np.any((distance!=0)&(np.abs(distance)<minimum_vertex_distance_m)):
        raise ValueError('Choose a plane outside the declared vertex separation guard')
    d=distance[f];cross=(d.min(1)<0)&(d.max(1)>0)
    identity=np.eye(3);untouched=np.flatnonzero(~cross);generic=np.flatnonzero(cross & ~np.any(d==0,axis=1))
    through=np.flatnonzero(cross & np.any(d==0,axis=1))
    if not len(generic) and not len(through):
        return dict(positions=p.copy(),faces=f.copy(),parent_faces=np.arange(len(f)),
            corner_barycentrics=np.broadcast_to(identity,(len(f),3,3)).copy(),
            evidence=dict(crossed_source_faces=0,new_controls=0,added_faces=0,source_surface_unchanged=True))
    # Rotate each source face; a cyclic permutation preserves its orientation.
    dg=d[generic];lone=np.where(np.sum(dg>0,axis=1)==1,np.argmax(dg,axis=1),np.argmin(dg,axis=1))
    order=(lone[:,None]+np.arange(3))%3;fg=np.take_along_axis(f[generic],order,axis=1)
    values=np.take_along_axis(dg,order,axis=1)
    zero=np.argmax(d[through]==0,axis=1);zo=(zero[:,None]+np.arange(3))%3
    ft=np.take_along_axis(f[through],zo,axis=1);vd=np.take_along_axis(d[through],zo,axis=1)
    edges=np.r_[fg[:,[0,1]],fg[:,[0,2]],ft[:,[1,2]]]
    unique,edge_inverse=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True)
    a,b=unique.T;parameter=distance[a]/(distance[a]-distance[b])
    intersections=p[a]+parameter[:,None]*(p[b]-p[a]);intersections[:,axis]=plane
    vertices=np.r_[p,intersections];edge_inverse+=len(p);ng=len(generic)
    ab=edge_inverse[:ng];ac=edge_inverse[ng:2*ng];bc=edge_inverse[2*ng:]
    ga,gb,gc=fg.T
    gf=np.stack([np.c_[ga,ab,ac],np.c_[ab,gb,gc],np.c_[ab,gc,ac]],axis=1).reshape(-1,3)
    ba,bb,cc=identity[order[:,0]],identity[order[:,1]],identity[order[:,2]]
    tab=values[:,0]/(values[:,0]-values[:,1]);tac=values[:,0]/(values[:,0]-values[:,2])
    bab=ba+tab[:,None]*(bb-ba);bac=ba+tac[:,None]*(cc-ba)
    gbary=np.stack([np.stack([ba,bab,bac],axis=1),np.stack([bab,bb,cc],axis=1),
                   np.stack([bab,cc,bac],axis=1)],axis=1).reshape(-1,3,3)
    ta,tb,tc=ft.T
    tf=np.stack([np.c_[ta,tb,bc],np.c_[ta,bc,tc]],axis=1).reshape(-1,3)
    za,zb,zc=identity[zo[:,0]],identity[zo[:,1]],identity[zo[:,2]]
    tbc=vd[:,1]/(vd[:,1]-vd[:,2]);zbc=zb+tbc[:,None]*(zc-zb)
    tbary=np.stack([np.stack([za,zb,zbc],axis=1),np.stack([za,zbc,zc],axis=1)],axis=1).reshape(-1,3,3)
    result=np.r_[f[untouched],gf,tf]
    parents=np.r_[untouched,np.repeat(generic,3),np.repeat(through,2)]
    bary=np.concatenate([np.broadcast_to(identity,(len(untouched),3,3)),gbary,tbary])
    source=p[f[parents]];target=vertices[result]
    sn=np.cross(source[:,1]-source[:,0],source[:,2]-source[:,0])
    tn=np.cross(target[:,1]-target[:,0],target[:,2]-target[:,0])
    if np.any(np.einsum('ni,ni->n',sn,tn)<=0):
        raise ValueError('The proposed plane produced a collapsed or reversed physical facet')
    return dict(positions=vertices,faces=result,parent_faces=parents,corner_barycentrics=bary,
        evidence=dict(crossed_source_faces=int(cross.sum()),new_controls=len(unique),
            added_faces=len(result)-len(f),source_surface_unchanged=True,
            original_controls_literal=True,shared_edges_split_once=True,
            corner_attribute_ancestry=True,geometric_storage_and_anatomical_checks_pending=True))
