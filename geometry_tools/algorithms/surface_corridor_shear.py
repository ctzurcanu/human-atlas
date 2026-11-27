"""Fit a surface into a measured Z corridor without changing its XY chart.

Inputs use one consistent unit. A cubic sheet displacement z'=z+d(x,y) has
continuous Jacobian determinant one. Discrete exported triangles must still
be checked independently for crossings, volume and finite containment.
"""
import numpy as np
import scipy.sparse as sp
from scipy.linalg import qr
from scipy.optimize import linprog
import clarabel


def cross2(a,b):return a[0]*b[1]-a[1]*b[0]


def intersection(a,b):
    poly=[p.copy()for p in a];sign=np.sign(cross2(b[1]-b[0],b[2]-b[0]))
    for i in range(3):
        if not poly:break
        p,q=b[i],b[(i+1)%3];edge=q-p;out=[]
        for j,current in enumerate(poly):
            previous=poly[j-1];dc=sign*cross2(edge,current-p);dp=sign*cross2(edge,previous-p)
            if(dc>=0)!=(dp>=0):out.append(previous+dp/(dp-dc)*(current-previous))
            if dc>=0:out.append(current)
        poly=out
    return np.array(poly)


def bary(points,tri):
    uv=np.linalg.solve(np.column_stack([tri[1]-tri[0],tri[2]-tri[0]]),(points-tri[0]).T).T
    return np.column_stack([1-uv.sum(1),uv])


def cubic(d):
    r=abs(d)
    return np.where(r<1,(4-6*r*r+3*r**3)/6,np.where(r<2,(2-r)**3/6,0.))


def cubic_derivative(d):
    r=abs(d)
    return np.sign(d)*np.where(r<1,-2*r+1.5*r*r,np.where(r<2,-.5*(2-r)**2,0.))


def basis(points,knots,spacing):
    d=(points[:,:2,None].transpose(0,2,1)-knots[None])/spacing
    return cubic(d[:,:,0])*cubic(d[:,:,1])


def gradient(points,knots,coefficients,spacing):
    d=(points[:,:2,None].transpose(0,2,1)-knots[None])/spacing
    return np.column_stack([(cubic_derivative(d[:,:,0])*cubic(d[:,:,1]))@coefficients,
                            (cubic(d[:,:,0])*cubic_derivative(d[:,:,1]))@coefficients])/spacing


def fit_vertex_corridor(vertices,faces,behind,ahead,*,pinned=None,guard=0.,minimum_self_fraction=1.,moving_faces_only=False,preserve_peers=(),declared_contact_peers=()):
    """Constrain the actual exported linear triangles, keeping XY literal.

    Unique XY controls share a Z displacement. The finite triangular mesh is
    checked by the caller; this discretization does not claim an analytical
    volume-preserving spline. Exact projected polygon vertices suffice for
    linear triangle-pair inequalities over their entire overlap.

    ``moving_faces_only`` excludes only triangles whose three XY controls are
    explicitly pinned. Every triangle that can move still receives all peer
    constraints. This permits a local refit without claiming to repair old
    violations on unchanged parts of a reference surface; callers must retain
    those limitations and independently check all new contacts.

    ``preserve_peers`` keeps each originally separated projected triangle pair
    on its original side. It does not assume that an entire folded neighbour
    lies in front of or behind the source. Pre-existing crossings must first
    be repaired or pinned; ambiguous pairs are reported instead of accepted.

    A caller may explicitly name already adjacent ``declared_contact_peers``.
    Only zero-gap coplanar facets with opposed outward normals can then use
    the declared material side. Their gap remains nonnegative; this does not
    infer an attachment or permit crossing a touching neighbour. Same-facing
    and genuinely crossing pairs remain ambiguous.
    """
    v=np.asarray(vertices);f=np.asarray(faces);tri=v[f]
    declared_contacts=set(declared_contact_peers)
    if not declared_contacts.issubset({name for name,_ in preserve_peers}):
        raise ValueError('Declared contacts must name explicit preserved peers')
    xy,inv=np.unique(v[:,:2],axis=0,return_inverse=True)
    locked=np.zeros(len(xy),bool)
    if pinned is not None:locked[inv[pinned]]=True
    active=np.flatnonzero(~locked);lookup=np.full(len(xy),-1);lookup[active]=np.arange(len(active))
    movable=np.any(lookup[inv[f]]>=0,axis=1)
    rows=[];cols=[];values=[];bounds=[];counts={};nz=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])[:,2]
    for lower,peers in [(True,behind),(False,ahead)]:
        own=np.flatnonzero((nz<0 if lower else nz>0)&(movable if moving_faces_only else True))
        for name,other in peers:
            other=np.asarray(other);oz=np.cross(other[:,1]-other[:,0],other[:,2]-other[:,0])[:,2]
            other=other[oz>0 if lower else oz<0];pl,ph=other[:,:,:2].min(1),other[:,:,:2].max(1);count=0
            for fi in own:
                t=tri[fi];lo,hi=t[:,:2].min(0),t[:,:2].max(0)
                for oi in np.flatnonzero(np.all(ph>=lo,axis=1)&np.all(pl<=hi,axis=1)):
                    peer=other[oi];poly=intersection(t[:,:2],peer[:,:2])
                    if len(poly)<3:continue
                    weights=bary(poly,t[:,:2]);pw=bary(poly,peer[:,:2]);delta=pw@peer[:,2]-weights@t[:,2]
                    indices=lookup[inv[f[fi]]]
                    for wi,di in zip(weights,delta):
                        ri=len(bounds);bounds.append(float(-di-guard if lower else di-guard))
                        for column,value in zip(indices,wi):
                            if column>=0 and value!=0:
                                rows.append(ri);cols.append(int(column));values.append(float(-value if lower else value))
                    count+=len(poly)
            counts[name]=count
    for name,other in preserve_peers:
        other=np.asarray(other);oz=np.cross(other[:,1]-other[:,0],other[:,2]-other[:,0])[:,2]
        usable=np.flatnonzero(oz!=0);pl,ph=other[:,:,:2].min(1),other[:,:,:2].max(1);count=0;contact_count=0
        for fi in np.flatnonzero((nz!=0)&(movable if moving_faces_only else True)):
            t=tri[fi];lo,hi=t[:,:2].min(0),t[:,:2].max(0)
            eligible=usable[np.all(ph[usable]>=lo,axis=1)&np.all(pl[usable]<=hi,axis=1)]
            for oi in eligible:
                peer=other[oi];poly=intersection(t[:,:2],peer[:,:2])
                if len(poly)<3:continue
                area=sum(cross2(poly[i]-poly[0],poly[i+1]-poly[0])for i in range(1,len(poly)-1))
                if abs(area)<1e-16:continue
                mw,pw=bary(poly,t[:,:2]),bary(poly,peer[:,:2]);gap=mw@t[:,2]-pw@peer[:,2]
                eps=max(1.,float(abs(np.r_[t,peer]).max()))*np.finfo(float).eps*128
                if gap.min()>=-eps and gap.max()>eps:direction=1.
                elif gap.max()<=eps and gap.min()<-eps:direction=-1.
                elif name in declared_contacts and np.max(abs(gap))<=eps:
                    own_normal=np.cross(t[1]-t[0],t[2]-t[0]);peer_normal=np.cross(peer[1]-peer[0],peer[2]-peer[0])
                    cosine=float(own_normal@peer_normal/(np.linalg.norm(own_normal)*np.linalg.norm(peer_normal)))
                    if cosine>-1+1e-8:
                        return dict(status='OriginalPeerOrderAmbiguous',peer=name,own_face=int(fi),peer_face=int(oi),
                                    original_gap_range=[float(gap.min()),float(gap.max())],counts=counts,
                                    declared_contact_rejected='Coplanar material normals are not opposed')
                    direction=-float(np.sign(nz[fi]));contact_count+=1
                else:
                    return dict(status='OriginalPeerOrderAmbiguous',peer=name,own_face=int(fi),peer_face=int(oi),
                                original_gap_range=[float(gap.min()),float(gap.max())],counts=counts)
                indices=lookup[inv[f[fi]]]
                for wi,g in zip(mw,gap):
                    ordered_gap=max(0.,float(direction*g))
                    ri=len(bounds);bounds.append(ordered_gap-min(guard,ordered_gap))
                    for column,value in zip(indices,wi):
                        if column>=0 and value!=0:
                            rows.append(ri);cols.append(int(column));values.append(float(-direction*value))
                    count+=1
        counts[name+'-preserved-order']=count
        if contact_count:counts[name+'-declared-zero-gap-contacts']=contact_count
    # Preserve the order of every pair of projected source surface triangles.
    # A fraction below one explicitly permits local profile compression; it
    # is a construction bound, never an intersection tolerance or anatomical
    # measurement. Shared edges retain their zero gap exactly.
    assert 0<minimum_self_fraction<=1
    pl,ph=tri[:,:,:2].min(1),tri[:,:,:2].max(1);self_count=0
    for fi,t in enumerate(tri):
        if nz[fi]==0:continue
        for oi in np.flatnonzero(np.all(ph>=pl[fi],axis=1)&np.all(pl<=ph[fi],axis=1)&(np.arange(len(f))>fi)):
            peer=tri[oi]
            if nz[oi]==0:continue
            poly=intersection(t[:,:2],peer[:,:2])
            if len(poly)<3:continue
            area=sum(cross2(poly[i]-poly[0],poly[i+1]-poly[0])for i in range(1,len(poly)-1))
            if abs(area)<1e-16:continue
            mw,pw=bary(poly,t[:,:2]),bary(poly,peer[:,:2]);gap=pw@peer[:,2]-mw@t[:,2]
            scale=max(1.,float(abs(tri).max()));eps=scale*np.finfo(float).eps*128
            if gap.min()>=-eps and gap.max()>eps:direction=1.
            elif gap.max()<=eps and gap.min()<-eps:direction=-1.
            else:continue
            for own_w,peer_w,g in zip(mw,pw,gap):
                # A common source point can interpolate to a signed residual
                # of a few double ULPs. Its exact shared-coordinate gap is zero.
                ordered_gap=max(0.,float(direction*g))
                ri=len(bounds);bounds.append(float((1-minimum_self_fraction)*ordered_gap))
                for sign,weights,face in [(direction,own_w,f[fi]),(-direction,peer_w,f[oi])]:
                    for column,value in zip(lookup[inv[face]],weights):
                        if column>=0 and value!=0:rows.append(ri);cols.append(int(column));values.append(float(sign*value))
                self_count+=1
    counts['source-surface-order']=self_count
    A=sp.csc_matrix((values,(rows,cols)),shape=(len(bounds),len(active)));b=np.array(bounds)
    # Near a pinned corner, projected constraints can have vanishing response
    # and large positive slack. An explicit displacement box, derived from
    # the full input Z extent, makes those rows provably redundant and avoids
    # arbitrarily large normalized RHS values during solver equilibration.
    all_z=[v[:,2],*[np.asarray(t)[:,:,2].ravel()for _,t in behind+ahead]]
    maximum_movement=float(max(x.max()for x in all_z)-min(x.min()for x in all_z))
    assert maximum_movement>0
    redundant=b>=np.asarray(abs(A).sum(1)).ravel()*maximum_movement
    A=A[~redundant];b=b[~redundant]
    A=sp.vstack([A,sp.eye(len(active)),-sp.eye(len(active))],format='csc')
    b=np.r_[b,np.full(2*len(active),maximum_movement)]
    row_norm=np.sqrt(np.asarray(A.multiply(A).sum(1)).ravel());zero=row_norm==0
    if np.any(b[zero]<0):return dict(status='InfeasibleFixedTriangle',counts=counts,minimumFixedSlack=float(b[zero].min()))
    A=sp.diags(1/row_norm[~zero])@A[~zero];b=b[~zero]/row_norm[~zero]
    edges=np.unique(np.sort(inv[f[:,[[0,1],[1,2],[2,0]]]].reshape(-1,2),axis=1),axis=0)
    edges=edges[edges[:,0]!=edges[:,1]];rr=[];cc=[];dd=[]
    for ri,edge in enumerate(edges):
        for sign,vi in zip([-1.,1.],edge):
            ci=lookup[vi]
            if ci>=0:rr.append(ri);cc.append(int(ci));dd.append(sign)
    D=sp.csc_matrix((dd,(rr,cc)),shape=(len(edges),len(active)))
    P=sp.eye(len(active),format='csc')+.15*(D.T@D)
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.max_iter=100;settings.time_limit=15.
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-10
    P=sp.triu(P,format='csc');A=sp.csc_matrix(A)
    for matrix in [P,A]:matrix.sum_duplicates();matrix.sort_indices()
    original_b=b.copy();original_A=A.copy()
    result=clarabel.DefaultSolver(P,np.zeros(len(active)),A,b,
                                  [clarabel.NonnegativeConeT(len(b))],settings).solve()
    if str(result.status)!='Solved':
        feasible=linprog(np.zeros(len(active)),A_ub=A,b_ub=b,bounds=(None,None),method='highs',options={'time_limit':10.})
        return dict(status=str(result.status),counts=counts,unknowns=len(active),inequalities=len(b),
                    independentLinearFeasibility=bool(feasible.success),linearSolverMessage=feasible.message)
    d=np.zeros(len(xy));d[active]=result.x;moved=v.copy();moved[:,2]+=d[inv]
    residual=float(np.min(original_b-original_A@result.x))
    assert residual>=-1e-8,(str(result.status),residual,result.r_prim,result.r_dual)
    return dict(status='Solved',counts=counts,positions=moved,displacements=d[inv],
                minimumConstraintSlack=residual,uniqueXY=xy,displacementControls=d)


def fit_z_corridor(vertices,faces,behind,ahead,*,spacing,pinned=None,guard=0.,posterior_only=True,diagnostic_path=None,
                   knot_bounds=None,nonpositive_vertices=False,coefficient_bound=None,displacement_weights=None,solver_variable_scale=1.,moving_faces_only=False,peer_guards=None,triangle_pairs=None):
    """Use actual projected triangle polygons, not bounding-box clearance.

    behind/ahead contain (name, triangles) pairs. Behind surfaces constrain the
    posterior muscle envelope; ahead surfaces constrain its anterior envelope.
    The caller must establish these anatomical sides and handle multi-interval
    geometry. The routine never closes an opening or invents a tissue volume.

    Optional finite weights in [0,1] multiply the spline response at each
    source control. They permit a medial contour correction while the opposite
    surface stays fixed. Such a weighted field is not a unit-Jacobian shear:
    callers must provide its continuous extension, derivative and positive
    Jacobian bound, as well as checking the actual exported triangles.

    ``moving_faces_only`` omits only facets with literally zero spline
    response at all three controls. It is a local refit, not validation of
    unchanged inherited contacts; before/after whole-surface gates remain
    the caller's responsibility.

    Explicit ``triangle_pairs`` may restrict a named compound peer to measured
    contact pairs. This avoids assigning one global side to unrelated folds
    of that peer. It is only a local contact solve; callers must still check
    all emitted surfaces against the entire peer after applying the field.
    """
    v=np.asarray(vertices);tri=v[faces]
    peer_guards={}if peer_guards is None else dict(peer_guards)
    if set(peer_guards)-{name for name,_ in behind+ahead}or any(not np.isfinite(x)or x<0 for x in peer_guards.values()):
        raise ValueError('Peer guards must name actual peers and be finite nonnegative distances')
    triangle_pairs={}if triangle_pairs is None else dict(triangle_pairs)
    peers_by_name={name:np.asarray(t)for name,t in behind+ahead}
    if set(triangle_pairs)-set(peers_by_name):raise ValueError('Pair restrictions must name actual peers')
    restrictions={}
    for name,pairs in triangle_pairs.items():
        pairs=np.asarray(pairs)
        if pairs.ndim!=2 or pairs.shape[1]!=2 or pairs.dtype.kind not in 'iu'or np.any(pairs<0)or np.any(pairs[:,0]>=len(faces))or np.any(pairs[:,1]>=len(peers_by_name[name])):
            raise ValueError('Pair restrictions must identify actual source/peer triangle indices')
        restrictions[name]={int(i):set(map(int,pairs[pairs[:,0]==i,1]))for i in np.unique(pairs[:,0])}
    if not np.isfinite(solver_variable_scale)or solver_variable_scale<=0:raise ValueError('Positive finite solver variable scale required')
    if coefficient_bound is not None and(not np.isfinite(coefficient_bound)or coefficient_bound<=0):raise ValueError('Positive finite coefficient bound required')
    gx=np.arange(np.floor(v[:,0].min()/spacing)-2,np.ceil(v[:,0].max()/spacing)+3)*spacing
    gy=np.arange(np.floor(v[:,1].min()/spacing)-2,np.ceil(v[:,1].max()/spacing)+3)*spacing
    knots=np.array([(x,y)for x in gx for y in gy])
    if knot_bounds is not None:
        bound=np.asarray(knot_bounds,float)
        if bound.shape!=(2,2)or not np.isfinite(bound).all()or np.any(bound[0]>=bound[1]):raise ValueError('Finite increasing 2D knot bounds required')
        knots=knots[np.all(knots>=bound[0],axis=1)&np.all(knots<=bound[1],axis=1)]
        if not len(knots):raise ValueError('Declared knot support contains no controls')
    B=basis(v,knots,spacing)
    if displacement_weights is not None:
        weights=np.asarray(displacement_weights,float)
        if weights.shape!=(len(v),)or not np.isfinite(weights).all()or np.any((weights<0)|(weights>1)):
            raise ValueError('One finite displacement weight in [0,1] per source control is required')
        B*=weights[:,None]
    pinned=np.zeros(len(v),dtype=bool)if pinned is None else pinned
    active=~np.any(B[pinned]>0,axis=0)if posterior_only and np.any(pinned)else np.ones(len(knots),dtype=bool)
    assert np.any(active)
    B=B[:,active];knots=knots[active]
    normal_z=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])[:,2]
    movable=np.any(np.any(B[faces]!=0,axis=2),axis=1)
    matrices=[];bounds=[];counts={};row_peers=[];row_faces=[]
    for lower,peers in [(True,behind),(False,ahead)]:
        own=np.flatnonzero((normal_z<0 if lower else normal_z>0)&(movable if moving_faces_only else True))
        for name,other in peers:
            peer_guard=peer_guards.get(name,guard)
            other=np.asarray(other)
            oz=np.cross(other[:,1]-other[:,0],other[:,2]-other[:,0])[:,2]
            original_peer_ids=np.flatnonzero(oz>0 if lower else oz<0);other=other[original_peer_ids]
            pl,ph=other[:,:,:2].min(1),other[:,:,:2].max(1);count=0
            for fi in own:
                allowed=restrictions.get(name)
                if allowed is not None and int(fi)not in allowed:continue
                t=tri[fi];low,high=t[:,:2].min(0),t[:,:2].max(0)
                eligible=np.flatnonzero(np.all(ph>=low,axis=1)&np.all(pl<=high,axis=1))
                if allowed is not None:eligible=eligible[np.array([int(original_peer_ids[j])in allowed[int(fi)]for j in eligible],bool)]
                for oi in eligible:
                    peer=other[oi];poly=intersection(t[:,:2],peer[:,:2])
                    if len(poly)<3:continue
                    mw,pw=bary(poly,t[:,:2]),bary(poly,peer[:,:2])
                    response=mw@B[faces[fi]]
                    delta=pw@peer[:,2]-mw@t[:,2]
                    if lower:
                        matrices.append(-response);bounds.append(-delta-peer_guard)
                    else:
                        matrices.append(response);bounds.append(delta-peer_guard)
                    row_peers.extend([name]*len(poly));row_faces.extend([int(fi)]*len(poly))
                    count+=len(poly)
            counts[name]=count
    A=np.vstack(matrices);b=np.concatenate(bounds)
    # A quadratic movement/bending objective chooses the least local change.
    distance=np.linalg.norm(knots[:,None]-knots[None],axis=2)
    edges=np.argwhere(np.triu(abs(distance-spacing)<spacing*1e-8,1))
    rr=np.repeat(np.arange(len(edges)),2);cc=edges.ravel()
    D=sp.csc_matrix((np.tile([1.,-1.],len(edges)),(rr,cc)),shape=(len(edges),len(knots)))
    P=sp.csc_matrix(B.T@B)+.05*(D.T@D)+sp.eye(len(knots),format='csc')*1e-6
    if posterior_only:
        A=np.vstack([A,np.eye(len(knots))]);b=np.r_[b,np.zeros(len(knots))]
        row_peers.extend(['coefficient-direction']*len(knots));row_faces.extend([-1]*len(knots))
    if nonpositive_vertices:
        # Signed coefficients can satisfy a fixed point inside the support;
        # the actual supplied controls still move only towards the interior.
        A=np.vstack([A,B]);b=np.r_[b,np.zeros(len(v))]
        row_peers.extend(['control-direction']*len(v));row_faces.extend([-1]*len(v))
    all_z=[v[:,2],*[np.asarray(t)[:,:,2].ravel()for _,t in behind+ahead]]
    maximum_movement=float(max(x.max()for x in all_z)-min(x.min()for x in all_z))
    if coefficient_bound is not None:maximum_movement=min(maximum_movement,float(coefficient_bound))
    redundant=b>=abs(A).sum(1)*maximum_movement
    row_peers=np.asarray(row_peers)[~redundant];row_faces=np.asarray(row_faces)[~redundant]
    A=A[~redundant];b=b[~redundant]
    A=np.vstack([A,np.eye(len(knots)),-np.eye(len(knots))])
    b=np.r_[b,np.full(2*len(knots),maximum_movement)]
    row_peers=np.r_[row_peers,np.repeat('coefficient-box',2*len(knots))];row_faces=np.r_[row_faces,np.full(2*len(knots),-1)]
    # Equality at a pinned point does not require each nearby signed spline
    # coefficient to vanish. Keep the actual linear equality and eliminate
    # only numerically dependent rows using a rank-revealing factorization.
    equality=np.empty((0,len(knots)))
    if not posterior_only and np.any(pinned):
        raw=B[pinned]
        _,R,pivot=qr(raw.T,pivoting=True,mode='economic')
        diagonal=abs(np.diag(R))
        rank=int(np.sum(diagonal>max(raw.shape)*np.finfo(float).eps*diagonal.max()))
        equality=raw[pivot[:rank]]
    row_scale=np.linalg.norm(A,axis=1)
    zero=row_scale==0
    if np.any(b[zero]<0):
        ids=np.flatnonzero(zero&(b<0));ids=ids[np.argsort(b[ids])]
        return dict(status='InfeasibleFixedRow',counts=counts,knots=knots,coefficients=None,
                    fixedConstraintExamples=[dict(peer=str(row_peers[i]),owner_face=int(row_faces[i]),physical_slack=float(b[i]))for i in ids[:20]])
    A=A[~zero]/row_scale[~zero,None];b=b[~zero]/row_scale[~zero]
    if diagnostic_path is not None:
        np.savez_compressed(diagnostic_path,A=A,b=b,equality=equality,P=P.toarray(),
                            B=B,knots=knots,vertices=v,faces=faces,pinned=pinned,solver_variable_scale=np.array(solver_variable_scale))
    settings=clarabel.DefaultSettings();settings.verbose=False;settings.max_iter=300
    settings.time_limit=15.
    settings.tol_feas=settings.tol_gap_abs=settings.tol_gap_rel=1e-10
    # y=c/scale makes metre-sized coefficients order one. The objective is
    # multiplied by one positive constant only; physical constraints and
    # their final residual checks are unchanged.
    all_A=sp.csc_matrix(np.vstack([A,equality]));all_b=np.r_[b/solver_variable_scale,np.zeros(len(equality))]
    cones=[clarabel.NonnegativeConeT(len(b))]
    if len(equality):cones.append(clarabel.ZeroConeT(len(equality)))
    P=sp.triu(P,format='csc')
    for matrix in [P,all_A]:matrix.sum_duplicates();matrix.sort_indices()
    result=clarabel.DefaultSolver(P,np.zeros(len(knots)),all_A,all_b,cones,settings).solve()
    if str(result.status)!='Solved':
        feasibility=linprog(np.zeros(len(knots)),A_ub=A,b_ub=b,A_eq=equality if len(equality)else None,
                            b_eq=np.zeros(len(equality))if len(equality)else None,bounds=(None,None),method='highs',options={'time_limit':10.})
        return dict(status=str(result.status),counts=counts,knots=knots,coefficients=None,
                    independentLinearFeasibility=bool(feasibility.success),linearSolverMessage=feasibility.message,
                    unknowns=len(knots),inequalities=len(b),equalities=len(equality),solverVariableScale=solver_variable_scale,
                    solverPrimalResidual=float(result.r_prim),solverDualResidual=float(result.r_dual),
                    physicalMinimumConstraintSlack=float(np.min(b-A@(np.array(result.x)*solver_variable_scale))))
    coefficients=np.array(result.x)*solver_variable_scale;delta=B@coefficients
    assert float(np.min(b-A@coefficients))>=-1e-8
    delta[pinned]=0.
    moved=v.copy();moved[:,2]+=delta
    return dict(status=str(result.status),counts=counts,knots=knots,coefficients=coefficients,
                positions=moved,displacements=delta,gradient=gradient(v,knots,coefficients,spacing)if displacement_weights is None else None,
                minimumConstraintSlack=float(np.min(b-A@coefficients)),
                continuousJacobianDeterminant=1. if displacement_weights is None else None,
                weightedFieldDerivativeRequired=displacement_weights is not None,solverVariableScale=solver_variable_scale)
