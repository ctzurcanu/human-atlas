"""Choose faithful Float32 rounding without fitting or deforming anatomy.

Each mutable coordinate is restricted to its two bracketing Float32 values.
Source separating axes become linear integer constraints. Exact triangle
predicates, source orientation and duplicate-control checks reject proposals;
this helper does not replace a complete emitted-mesh or volume audit.
"""
import itertools,time
import numpy as np
from .tube_sections import TubeGeometryError
from ..distance import triangle_pairs


def oriented_rounding(source, encoded, faces, corner_words, *, origin=(0.,0.,0.),
                      mutable_controls=None, maximum_states=262144, time_limit=5.):
    """Repair diagnosed Float32 facet flips using only bracketing grid values.

    Small connected sets are enumerated exactly. Every incident facet and its
    three authored corner normals must remain outward. This is an encoding
    candidate; full surface, lumen and occupied-volume audits remain required.
    """
    p=np.asarray(source,float); q=np.array(encoded,float,copy=True)
    f=np.asarray(faces,np.int64); words=np.asarray(corner_words,float); o=np.asarray(origin,float)
    if p.shape!=q.shape or p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all()or not np.isfinite(q).all():
        raise TubeGeometryError('Corresponding finite source and encoded banks required')
    if f.ndim!=2 or f.shape[1]!=3 or words.shape!=(len(f),3,3):
        raise TubeGeometryError('Indexed triangles and their complete corner normal fields required')
    source_tri=p[f]; normal=np.cross(source_tri[:,1]-source_tri[:,0],source_tri[:,2]-source_tri[:,0])
    if not np.all(np.einsum('fci,fi->fc',words,normal)>0):
        raise TubeGeometryError('Repair the source orientation before Float32 encoding')
    tri=q[f]; encoded_normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0])
    bad=np.flatnonzero(np.any(np.einsum('fci,fi->fc',words,encoded_normal)<=0,axis=1))
    start=time.monotonic(); initial=q.copy()
    if not len(bad):
        return q,dict(faithful_floor_ceil_only=True,changed_controls=[],diagnosed_faces=0,
                      component_states=[],elapsed_seconds=time.monotonic()-start,full_mesh_acceptance=False)
    allowed=set(np.unique(f[bad]).tolist())if mutable_controls is None else set(map(int,mutable_controls))
    nearest=(p-o).astype(np.float32); decoded=nearest.astype(float)+o
    lower=np.where(decoded>p,np.nextafter(nearest,np.float32(-np.inf)),nearest).astype(float)+o
    upper=np.where(decoded<p,np.nextafter(nearest,np.float32(np.inf)),nearest).astype(float)+o
    selected=set(map(int,np.unique(f[bad])))&allowed
    groups=[]
    while selected:
        group={selected.pop()}
        while True:
            incident=np.isin(f,list(group)).any(1)
            more=set(map(int,np.unique(f[incident])))&selected
            if not more:break
            selected-=more;group|=more
        groups.append(sorted(group))
    states=[]; incident_counts=[]
    for group in groups:
        if time.monotonic()-start>=time_limit:raise TimeoutError('Float32 orientation enumeration reached its deadline')
        variables=[(i,k)for i in group for k in range(3)if upper[i,k]>lower[i,k]]
        count=1<<len(variables)
        if not variables or count>maximum_states:
            raise TubeGeometryError('The faithful orientation component exceeds its finite enumeration budget')
        incident=np.flatnonzero(np.isin(f,group).any(1)); local=f[incident]; ids=np.unique(local)
        index={int(v):i for i,v in enumerate(ids)}; mapped=np.array([[index[int(i)]for i in face]for face in local])
        choices=np.broadcast_to(q[ids],(count,len(ids),3)).copy()
        bits=(np.arange(count,dtype=np.uint64)[:,None]>>np.arange(len(variables),dtype=np.uint64))&1
        for j,(i,k)in enumerate(variables):choices[:,index[i],k]=lower[i,k]+bits[:,j]*(upper[i,k]-lower[i,k])
        triangles=choices[:,mapped]; n=np.cross(triangles[:,:,1]-triangles[:,:,0],triangles[:,:,2]-triangles[:,:,0])
        ok=np.all(np.einsum('fci,sfi->sfc',words[incident],n)>0,axis=(1,2))
        ok&=np.all(np.einsum('sfi,fi->sf',n,normal[incident])>0,axis=1)
        if not np.any(ok):raise TubeGeometryError('No faithful rounding satisfies all incident facet orientations')
        cost=np.sum((choices-p[ids])**2,axis=(1,2));cost[~ok]=np.inf
        q[ids]=choices[int(np.argmin(cost))];states.append(count);incident_counts.append(len(incident))
    actual=q[f]; n=np.cross(actual[:,1]-actual[:,0],actual[:,2]-actual[:,0])
    if not np.all(np.einsum('fci,fi->fc',words,n)>0):raise TubeGeometryError('A Float32 corner remains inward')
    if not np.array_equal((q-o).astype(np.float32).astype(float)+o,q):raise TubeGeometryError('A coordinate is not actually representable')
    changed=np.flatnonzero(np.any(q!=initial,axis=1))
    if not set(map(int,changed))<=allowed:raise TubeGeometryError('An immutable source control changed')
    if np.any((q[changed]<lower[changed])|(q[changed]>upper[changed])):raise TubeGeometryError('A coordinate left its bracketing values')
    return q,dict(faithful_floor_ceil_only=True,changed_controls=changed.tolist(),diagnosed_faces=len(bad),
                  component_states=states,incident_faces_checked=incident_counts,
                  maximum_encoded_control_change_m=float(np.linalg.norm(q-initial,axis=1).max()),
                  maximum_control_error_m=float(np.linalg.norm(q-p,axis=1).max()),
                  elapsed_seconds=time.monotonic()-start,full_mesh_acceptance=False)


def separating_axis_options(a,b,limit=4):
    reference=a[0];a=a-reference;b=b-reference
    ea=np.roll(a,-1,axis=0)-a;eb=np.roll(b,-1,axis=0)-b
    na=np.cross(ea[0],ea[1]);nb=np.cross(eb[0],eb[1])
    axes=[na,nb]+[np.cross(x,y) for x in ea for y in eb]
    axes += [np.cross(n,e) for n in [na,nb] for e in [*ea,*eb]]
    choices=[]
    for axis in axes:
        length=np.linalg.norm(axis)
        if length==0:continue
        axis=axis/length;aa=a@axis;bb=b@axis
        forward=float(bb.min()-aa.max());reverse=float(aa.min()-bb.max())
        choices.append((forward,axis));choices.append((reverse,-axis))
    best_gap=max(row[0]for row in choices)
    scale=max(float(np.ptp(np.r_[a,b],axis=0).max()),np.finfo(float).tiny)
    roundoff=512*np.finfo(float).eps*scale
    # At a contact several axes have the same zero gap. A common face normal
    # constrains no in-plane overlap and can force an unnecessarily exact
    # coplanar grid embedding. Prefer an actual separating edge/vertex feature:
    # fewer controls lie on its contact plane, leaving a stable interior margin.
    def feature_score(row):
        gap,axis=row;aa=a@axis;bb=b@axis;plane=(aa.max()+bb.min())/2
        margins=np.r_[plane-aa,bb-plane];contact=int(np.sum(np.abs(margins)<=roundoff))
        positive=margins[margins>roundoff]
        interior=float(positive.min())if len(positive)else 0.
        return -contact,interior,gap
    if best_gap< -1e-16:raise TubeGeometryError('No stable source separating axis; repair the anatomical source first')
    tied=[row for row in choices if row[0]>=best_gap-roundoff]
    # At a vertex contact, sums of valid support directions enter the support
    # cone interior and avoid overconstraining unrelated coplanar controls.
    for first,second in itertools.combinations(tied,2):
        axis=first[1]+second[1];length=np.linalg.norm(axis)
        if length<1e-8:continue
        axis/=length;gap=float((b@axis).min()-(a@axis).max())
        if gap>=best_gap-roundoff:tied.append((gap,axis))
    result=[]
    for gap,axis in sorted(tied,key=feature_score,reverse=True):
        if feature_score((gap,axis))[0]==-6:continue
        if any(np.dot(axis,old[0])>1-1e-12 for old in result):continue
        result.append((axis,max(0.,gap)))
        if len(result)>=limit:break
    if not result:
        gap,axis=max(tied,key=feature_score);result=[(axis,max(0.,gap))]
    return result


def _separating_axis(a,b):
    return separating_axis_options(a,b,limit=1)[0]


def directed_rounding(source,encoded,faces,pairs,*,origin=(0,0,0),mutable_controls=None,time_limit=10.,maximum_rounds=5,corner_words=None,additional_controls=(),validation_pairs=None):
    """Return a candidate grid embedding, leaving ``source`` literally fixed.

    ``pairs`` is [N,2,3] in one global control bank. The caller supplies the
    actual diagnosed triangle pairs, never inferred anatomical connections.
    ``validation_pairs`` optionally retains an explicitly enumerated nearby
    triangle neighborhood using exact rejection, without forcing its source
    support planes onto already rounded, immutable neighboring controls.
    Remaining/new defects must be found by independent whole-mesh checks.
    """
    from scipy.optimize import milp,Bounds,LinearConstraint
    from scipy.sparse import coo_matrix
    orientation_only=pairs is None
    p=np.asarray(source,float);q=np.array(encoded,float,copy=True);f=np.asarray(faces,np.int64);pairs=np.empty((0,2,3),np.int64)if orientation_only else np.asarray(pairs,np.int64);o=np.asarray(origin,float)
    words=np.asarray(corner_words,float)if corner_words is not None else None
    if p.shape!=q.shape or p.ndim!=2 or p.shape[1]!=3 or not np.isfinite(p).all() or not np.isfinite(q).all():raise TubeGeometryError('Corresponding finite source and encoded banks required')
    if pairs.ndim!=3 or pairs.shape[1:]!=(2,3) or(not len(pairs)and not orientation_only):raise TubeGeometryError('Supply actual indexed triangle pairs')
    if orientation_only and(words is None or not len(additional_controls)):raise TubeGeometryError('Orientation-only rounding requires corner fields and explicitly diagnosed controls')
    if words is not None and (words.shape!=(len(f),3,3)or not np.isfinite(words).all()):raise TubeGeometryError('Complete finite corner normal fields required')
    retained=np.empty((0,2,3),np.int64)if validation_pairs is None else np.asarray(validation_pairs,np.int64)
    if retained.ndim!=3 or retained.shape[1:]!=(2,3):raise TubeGeometryError('Validation pairs must be explicit indexed triangle pairs')
    if retained.size and(retained.min()<0 or retained.max()>=len(p)):raise TubeGeometryError('A validation control is outside the source bank')
    checked=np.concatenate([pairs,retained])
    source_copy=p.copy();start=time.monotonic();relation=triangle_pairs(p[checked[:,0]],p[checked[:,1]])['relation']
    if np.isin(relation,['crossing','coplanar-overlap']).any():raise TubeGeometryError('The anatomical source already intersects')
    selected=np.unique(np.r_[pairs.ravel(),np.asarray(additional_controls,np.int64).ravel()])
    if len(selected)and(selected.min()<0 or selected.max()>=len(p)):raise TubeGeometryError('A selected control is outside the source bank')
    allowed=set(selected.tolist()) if mutable_controls is None else set(map(int,mutable_controls));selected=np.array([i for i in selected if int(i) in allowed],np.int64)
    nearest=(p-o).astype(np.float32);decoded=nearest.astype(float)+o
    lower=np.where(decoded>p,np.nextafter(nearest,np.float32(-np.inf)),nearest).astype(float)+o
    upper=np.where(decoded<p,np.nextafter(nearest,np.float32(np.inf)),nearest).astype(float)+o
    variables=[(int(i),k) for i in selected for k in range(3) if upper[i,k]>lower[i,k]]
    if not variables:raise TubeGeometryError('No mutable faithful rounding choices are available')
    columns={v:i for i,v in enumerate(variables)};base=q.copy();step=np.zeros_like(q)
    for i,k in variables:base[i,k]=lower[i,k];step[i,k]=upper[i,k]-lower[i,k]
    max_step=float(step.max());cost=np.array([((upper[i,k]-p[i,k])**2-(lower[i,k]-p[i,k])**2)/max_step**2 for i,k in variables])
    rows=[];bounds=[]
    def add(coeff,minimum=-np.inf,maximum=np.inf):
        coeff={k:v for k,v in coeff.items() if v!=0}
        if not coeff:
            if minimum>1e-12 or maximum< -1e-12:raise TubeGeometryError('Protected controls cannot satisfy the source separation')
            return
        scale=max(abs(v) for v in coeff.values());rows.append({k:v/scale for k,v in coeff.items()});bounds.append((minimum/scale,maximum/scale))
    for a,b in pairs:
        axis,gap=_separating_axis(p[a],p[b]);minimum=gap*.25
        for i,j in itertools.product(a,b):
            coeff={}
            for control,sign in [(int(i),-1),(int(j),1)]:
                for k in range(3):
                    key=columns.get((control,k))
                    if key is not None:coeff[key]=coeff.get(key,0.)+sign*axis[k]*step[control,k]
            add(coeff,minimum-float((base[j]-base[i])@axis))
    affected=np.flatnonzero(np.isin(f,selected).any(1));old=p[f[affected]];normal=np.cross(old[:,1]-old[:,0],old[:,2]-old[:,0]);length=np.linalg.norm(normal,axis=1)
    if np.any(length==0):raise TubeGeometryError('A source facet is degenerate')
    normal/=length[:,None];t=q[f[affected]];area=np.einsum('ij,ij->i',np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),normal)
    for ordinal,(face,n,vertices,source_area)in enumerate(zip(f[affected],normal,old,length)):
        constraints=[(n,source_area)]
        if words is not None:
            for word in words[affected[ordinal]]:
                size=np.linalg.norm(word)
                if size==0:raise TubeGeometryError('A source corner field is zero')
                direction=word/size;projection=source_area*float(n@direction)
                if projection<=0:raise TubeGeometryError('A source corner field is inward')
                constraints.append((direction,projection))
        for direction,minimum_area in constraints:
            gradient=np.array([np.cross(vertices[1]-vertices[2],direction),np.cross(vertices[2]-vertices[0],direction),np.cross(direction,vertices[1]-vertices[0])])
            coeff={};constant=float(np.cross(vertices[1]-vertices[0],vertices[2]-vertices[0])@direction)
            for corner,control in enumerate(face):
                for k in range(3):
                    column=columns.get((int(control),k))
                    if column is not None:
                        coeff[column]=coeff.get(column,0.)+gradient[corner,k]*step[control,k]
                        constant+=gradient[corner,k]*(base[control,k]-p[control,k])
            add(coeff,minimum_area*1e-4-constant)
    active=np.unique(f)
    def forbid(controls,bits):
        ids=sorted({columns[(int(i),k)] for i in controls for k in range(3) if (int(i),k) in columns})
        if not ids:raise TubeGeometryError('A rejected grid embedding has no mutable incident controls')
        add({j:1. if bits[j] else -1. for j in ids},maximum=float(sum(bits[j] for j in ids)-1))
    for attempt in range(maximum_rounds):
        rr=[];cc=[];data=[]
        for i,row in enumerate(rows):
            for j,value in row.items():rr.append(i);cc.append(j);data.append(value)
        matrix=coo_matrix((data,(rr,cc)),shape=(len(rows),len(variables))).tocsc()
        remaining=float(time_limit)-(time.monotonic()-start)
        if remaining<=0:raise TimeoutError('Directed Float32 solve reached its wall-clock deadline')
        result=milp(cost,integrality=np.ones(len(variables)),bounds=Bounds(np.zeros(len(variables)),np.ones(len(variables))),constraints=LinearConstraint(matrix,np.array([x[0] for x in bounds]),np.array([x[1] for x in bounds])),options={'time_limit':remaining,'mip_rel_gap':1e-9})
        if not result.success:
            error=TimeoutError if result.status==1 else TubeGeometryError
            raise error('Directed Float32 solve did not complete: '+str(result.message))
        bits=np.rint(result.x).astype(int);candidate=base.copy()
        for (i,k),bit in zip(variables,bits):candidate[i,k]+=step[i,k]*bit
        if not np.array_equal((candidate-o).astype(np.float32).astype(float)+o,candidate):raise TubeGeometryError('A proposed coordinate is not actually representable')
        actual=candidate[f[affected]];actual_normals=np.cross(actual[:,1]-actual[:,0],actual[:,2]-actual[:,0]);dot=np.einsum('ij,ij->i',actual_normals,normal);bad=np.flatnonzero(dot<=0)
        if words is not None:bad=np.union1d(bad,np.flatnonzero(np.any(np.einsum('fci,fi->fc',words[affected],actual_normals)<=0,axis=1)))
        _,inverse,count=np.unique(candidate[active],axis=0,return_inverse=True,return_counts=True);duplicates=np.flatnonzero(count>1)
        rechecked=triangle_pairs(candidate[checked[:,0]],candidate[checked[:,1]])['relation'];wrong=np.flatnonzero(np.isin(rechecked,['crossing','coplanar-overlap']))
        if not len(bad) and not len(duplicates) and not len(wrong):
            changed=np.flatnonzero(np.any(candidate!=q,axis=1));assert np.array_equal(p,source_copy)
            return candidate,dict(source_positions_literal=True,faithful_floor_ceil_only=True,orientation_only=orientation_only,changed_controls=changed.tolist(),selected_controls=selected.tolist(),variables=len(variables),constraints=len(rows),solver_rounds=attempt+1,solver_status=str(result.message),known_pairs_checked=len(pairs),validation_pairs_checked=len(retained),known_crossings_remaining=0,incident_faces_orientation_checked=len(affected),maximum_control_error_m=float(np.linalg.norm(candidate-p,axis=1).max()),maximum_encoded_control_change_m=float(np.linalg.norm(candidate-q,axis=1).max()),elapsed_seconds=time.monotonic()-start,full_mesh_acceptance=False)
        for i in bad:forbid(f[affected[i]],bits)
        for i in duplicates:forbid(active[inverse==i],bits)
        for i in wrong:forbid(checked[i].ravel(),bits)
    raise TubeGeometryError('Directed Float32 refinement budget exhausted')
