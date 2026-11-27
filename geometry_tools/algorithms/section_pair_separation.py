"""Sufficient triangle separation for shared rigid-section translation variables.

A caller supplies the route axis explicitly. This is a geometric halfspace,
not a rule for choosing anatomical sides or classifying a vessel junction.
"""
import numpy as np
from scipy import sparse


def directed_triangle_exit_axes(first,second,direction,*,clearance_m=0.):
    """SAT escape planes with least common translation along an explicit ray.

    For each literal triangle pair, chooses a full-interval separating plane
    that can be reached in the positive supplied direction. Subsequent variable
    section transport must enforce all nine vertex inequalities, not assume the
    triangles move rigidly together. This is not a whole-mesh collision proof.
    """
    a,b=np.asarray(first,float),np.asarray(second,float);direction=np.asarray(direction,float)
    if a.ndim!=3 or a.shape[1:]!=(3,3)or b.shape!=a.shape or direction.shape!=(3,)or not np.isfinite(a).all()or not np.isfinite(b).all()or not np.isfinite(direction).all()or np.linalg.norm(direction)<1e-12 or not np.isfinite(clearance_m)or clearance_m<0:
        raise ValueError('Finite matching triangles and an explicit positive escape ray required')
    direction=direction/np.linalg.norm(direction)
    ea=np.roll(a,-1,axis=1)-a;eb=np.roll(b,-1,axis=1)-b;na=np.cross(ea[:,0],ea[:,1]);nb=np.cross(eb[:,0],eb[:,1])
    axes=np.concatenate([na[:,None],nb[:,None],np.cross(ea[:,:,None],eb[:,None]).reshape(-1,9,3),np.cross(na[:,None],ea),np.cross(nb[:,None],eb)],axis=1)
    lengths=np.linalg.norm(axes,axis=2);axes=np.divide(axes,lengths[:,:,None],out=np.zeros_like(axes),where=lengths[:,:,None]>1e-25)
    rate=axes@direction;axes*=np.where(rate<0,-1.,1.)[:,:,None];rate=abs(rate)
    pa=np.einsum('pvc,pac->pav',a,axes);pb=np.einsum('pvc,pac->pav',b,axes);gap=pa.min(2)-pb.max(2)
    amount=np.divide(np.maximum(clearance_m-gap,0.),rate,out=np.full_like(rate,np.inf),where=rate>1e-10)
    chosen=amount.argmin(1);index=np.arange(len(a))
    if not np.isfinite(amount[index,chosen]).all():raise ValueError('No finite escape plane along the declared ray')
    return axes[index,chosen],amount[index,chosen]


def profile_radial_axes(first,second,first_stations,second_stations,pairs,*,preferred_axis):
    """Propose radial escape axes from literal source profile planes.

    Uses complete profile controls, rather than individual triangle normals.
    Protected non-profile core controls use adjacent source-station centres.
    Crossed profile tangents give the normal of their common course plane;
    parallel tangents use their radial centre offset. The declared preferred
    direction orients every proposal consistently across the whole route;
    individual face-centre differences must not flip the escape side. No anatomical side or
    collision-free path is established by this proposal.
    """
    pairs=np.asarray(pairs);preferred=np.asarray(preferred_axis,float)
    if pairs.ndim!=2 or pairs.shape[1]!=2 or preferred.shape!=(3,)or not np.isfinite(preferred).all()or np.linalg.norm(preferred)<1e-12:
        raise ValueError('Explicit source pairs and nonzero preferred direction required')
    def fields(mesh,station):
        station=np.asarray(station)
        if station.shape!=(len(mesh.vertices),)or station.dtype.kind not in 'iu'or np.any(station<0):raise ValueError('Original nonnegative profile ownership required')
        centre=np.zeros_like(mesh.vertices);tangent=np.zeros_like(mesh.vertices)
        groups={g:np.flatnonzero(station==g)for g in np.unique(station)}
        centres={g:mesh.vertices[ids].mean(0)for g,ids in groups.items()}
        for group,ids in groups.items():
            points=mesh.vertices[ids];c=centres[group]
            _,singular,v=np.linalg.svd(points-c,full_matrices=False)
            if len(ids)>=3 and singular[1]>1e-12:axis=v[-1]
            else:
                adjacent=np.unique(station[mesh.faces[np.any(np.isin(mesh.faces,ids),axis=1)]])
                offsets=np.array([centres[g]-c for g in adjacent if g!=group])
                if not len(offsets):raise ValueError('Unresolved source core has no neighbouring course')
                _,_,course=np.linalg.svd(offsets,full_matrices=False);axis=course[0]
            centre[ids]=c;tangent[ids]=axis
        return centre,tangent
    ca,ta=fields(first,first_stations);cb,tb=fields(second,second_stations)
    fa,fb=first.faces[pairs[:,0]],second.faces[pairs[:,1]]
    delta=ca[fa].mean(1)-cb[fb].mean(1)
    # Use the dominant literal profile on a face to avoid eigenvector sign
    # cancellation when a side triangle straddles adjacent rings.
    u=ta[fa[:,0]];v=tb[fb[:,0]];v*=np.where(np.einsum('ij,ij->i',u,v)<0,-1.,1.)[:,None]
    cross=np.cross(u,v);length=np.linalg.norm(cross,axis=1)
    aligned=u+v;aligned/=np.maximum(np.linalg.norm(aligned,axis=1),1e-20)[:,None]
    radial=delta-np.einsum('ij,ij->i',delta,aligned)[:,None]*aligned
    directions=np.where((length>.1)[:,None],cross,radial)
    lengths=np.linalg.norm(directions,axis=1);fallback=preferred-np.einsum('j,ij->i',preferred,aligned)[:,None]*aligned
    directions=np.where((lengths>1e-12)[:,None],directions,fallback)
    lengths=np.linalg.norm(directions,axis=1)
    if np.any(lengths<1e-12):raise ValueError('Preferred radial escape is parallel to an unresolved course')
    directions/=lengths[:,None]
    sign=directions@preferred
    ambiguous=abs(sign)<1e-8
    if ambiguous.any():
        replacement=fallback[ambiguous]
        magnitude=np.linalg.norm(replacement,axis=1)
        if np.any(magnitude<1e-12):raise ValueError('Preferred escape does not resolve the profile plane')
        directions[ambiguous]=replacement/magnitude[:,None]
        sign[ambiguous]=directions[ambiguous]@preferred
    directions*=np.where(sign<0,-1.,1.)[:,None]
    return directions


def triangle_axis_rows(first, second, first_weights, second_weights, pairs, axis, *, gap_m):
    w, z = np.asarray(first_weights, float), np.asarray(second_weights, float)
    pairs = np.asarray(pairs)
    axis = np.asarray(axis, float)
    if (w.ndim != 2 or z.ndim != 2 or w.shape[0] != len(first.vertices)
        or z.shape[0] != len(second.vertices) or w.shape[1] != z.shape[1]
        or not w.shape[1] or pairs.ndim != 2 or pairs.shape[1] != 2
        or pairs.dtype.kind not in 'iu' or not len(pairs)
        or axis.shape not in [(3,), (len(pairs), 3)] or not np.isfinite(axis).all()
        or not np.isfinite(w).all() or not np.isfinite(z).all()
        or not np.isfinite(gap_m) or gap_m < 0 or np.any(np.linalg.norm(axis,axis=-1) < 1e-12)
        or pairs.min() < 0 or pairs[:, 0].max() >= len(first.faces)
        or pairs[:, 1].max() >= len(second.faces) or first.frame != second.frame):
        raise ValueError('Explicit common-frame faces, translation weights, route axis and nonnegative gap required')
    directions = np.broadcast_to(axis, (len(pairs),3)).copy()
    directions /= np.linalg.norm(directions,axis=1)[:,None]
    fa, fb = first.faces[pairs[:, 0]], second.faces[pairs[:, 1]]
    weights = (w[fa][:, :, None, :] - z[fb][:, None, :, :]).reshape(-1, w.shape[1])
    base = (first.vertices[fa][:, :, None, :] - second.vertices[fb][:, None, :, :]).reshape(-1, 3)
    normal = np.repeat(directions,9,axis=0)
    rhs = gap_m - np.einsum('ij,ij->i',base,normal)
    return weights, normal, rhs


def separated_triangle_axes(first,second,*,prefer_coordinates=False):
    """Find the widest existing SAT gap, including coplanar edge normals.

    Floating projections propose conservative preservation planes. Positive
    margins must be replayed on actual emitted triangles; this is not an exact
    intersection predicate or a choice of anatomical routing for crossed pairs.
    """
    a,b=np.asarray(first,float),np.asarray(second,float)
    if a.ndim!=3 or a.shape[1:]!=(3,3)or b.shape!=a.shape or not np.isfinite(a).all()or not np.isfinite(b).all():
        raise ValueError('Matching finite literal triangle arrays required')
    if prefer_coordinates:
        positive=a.min(1)-b.max(1);negative=b.min(1)-a.max(1);gap=np.maximum(positive,negative)
        chosen=gap.argmax(1);index=np.arange(len(a));distance=gap[index,chosen]
        direction=np.eye(3)[chosen]*np.where(negative[index,chosen]>positive[index,chosen],-1.,1.)[:,None]
        unresolved=distance<=1e-9
        if unresolved.any():direction[unresolved],distance[unresolved]=separated_triangle_axes(a[unresolved],b[unresolved])
        return direction,distance
    ea=np.roll(a,-1,axis=1)-a;eb=np.roll(b,-1,axis=1)-b
    na=np.cross(ea[:,0],ea[:,1]);nb=np.cross(eb[:,0],eb[:,1])
    axes=np.concatenate([na[:,None],nb[:,None],np.cross(ea[:,:,None],eb[:,None]).reshape(-1,9,3),
        np.cross(na[:,None],ea),np.cross(nb[:,None],eb)],axis=1)
    length=np.linalg.norm(axes,axis=2);valid=length>1e-25
    axes=np.divide(axes,length[:,:,None],out=np.zeros_like(axes),where=valid[:,:,None])
    pa=np.einsum('pvc,pac->pav',a,axes);pb=np.einsum('pvc,pac->pav',b,axes)
    positive=pa.min(2)-pb.max(2);negative=pb.min(2)-pa.max(2)
    reverse=negative>positive;gap=np.maximum(positive,negative);gap[~valid]=-np.inf
    chosen=gap.argmax(1);index=np.arange(len(a));direction=axes[index,chosen]*np.where(reverse[index,chosen],-1.,1.)[:,None]
    return direction,gap[index,chosen]


def triangle_axis_sparse_rows(first,second,first_stations,second_stations,pairs,axis,*,station_count,gap_m):
    """Sparse equivalent for rigid station IDs; -1 denotes a fixed control.

    Identical station/axis rows retain only their strongest bound. This is an
    exact halfspace reduction, not geometric simplification or face removal.
    """
    a,b=np.asarray(first_stations),np.asarray(second_stations);pairs=np.asarray(pairs);axis=np.asarray(axis,float)
    gap=np.broadcast_to(np.asarray(gap_m,float),(len(pairs),))
    if (a.shape!=(len(first.vertices),)or b.shape!=(len(second.vertices),)or a.dtype.kind not in 'iu'or b.dtype.kind not in 'iu'
        or not isinstance(station_count,int)or station_count<1 or np.any(a < -1)or np.any(b < -1)
        or np.any(a>=station_count)or np.any(b>=station_count)or first.frame!=second.frame
        or pairs.ndim!=2 or pairs.shape[1]!=2 or pairs.dtype.kind not in 'iu'or not len(pairs)
        or pairs.min()<0 or pairs[:,0].max()>=len(first.faces)or pairs[:,1].max()>=len(second.faces)
        or axis.shape not in [(3,),(len(pairs),3)]or not np.isfinite(axis).all()
        or np.any(np.linalg.norm(axis,axis=-1)<1e-12)or not np.isfinite(gap).all()or np.any(gap<0)):
        raise ValueError('Explicit rigid stations, source pairs, route axes and finite nonnegative gaps required')
    fa,fb=first.faces[pairs[:,0]],second.faces[pairs[:,1]]
    direction=np.broadcast_to(axis,(len(pairs),3)).copy();direction/=np.linalg.norm(direction,axis=1)[:,None]
    normal=np.repeat(direction,9,axis=0)
    va=np.broadcast_to(fa[:,:,None],(len(pairs),3,3)).ravel();vb=np.broadcast_to(fb[:,None,:],(len(pairs),3,3)).ravel()
    rhs=np.repeat(gap,9)-np.einsum('ij,ij->i',first.vertices[va]-second.vertices[vb],normal)
    keys,inverse=np.unique(np.c_[a[va],b[vb],normal],axis=0,return_inverse=True)
    bounds=np.full(len(keys),-np.inf);np.maximum.at(bounds,inverse,rhs)
    rr=[];cc=[];vv=[]
    for column,sign in [(0,1.),(1,-1.)]:
        station=keys[:,column].astype(int)
        for component in range(3):
            use=(station>=0)&(keys[:,2+component]!=0)
            rr.append(np.flatnonzero(use));cc.append(station[use]*3+component);vv.append(sign*keys[use,2+component])
    matrix=sparse.coo_matrix((np.concatenate(vv),(np.concatenate(rr),np.concatenate(cc))),shape=(len(keys),3*station_count)).tocsr()
    matrix.eliminate_zeros()
    return matrix,bounds,dict(sourceTrianglePairs=len(pairs),inputVertexHalfspaces=len(rhs),strongestDistinctStationHalfspaces=len(bounds),
        nonzeroCoefficients=matrix.nnz,denseMatrixNotAllocated=True,sourceGeometryChanged=False)
