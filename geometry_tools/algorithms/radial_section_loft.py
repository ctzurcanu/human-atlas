"""Measured radial sections and coupled harmonic cross-section transport.

The caller supplies material curves and anatomy. A radial section is rejected
if it has multiple positive crossings; no hidden cavity is filled. Harmonic
interpolation alone is not a bijection certificate: every emitted section and
calculation tetrahedron must pass its orientation gate.
"""
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import splu


class PrismJacobianError(ValueError):
    def __init__(self,bad_slabs):
        super().__init__('A prism interpolation has a negative reference Jacobian')
        self.bad_slabs=list(map(int,bad_slabs))


def bound_section_deformations(reference,desired,faces,*,maximum_displacement_m,minimum_relative_area=.2):
    """Keep a shared section inside explicit displacement/Jacobian bounds.

    The reference is an affine placement of the same connected graph. This is
    a measured candidate fit, not permission to change an anatomical border.
    Every accepted triangle and the complete emitted volume still need checks.
    """
    baseline=np.asarray(reference,float);target=np.asarray(desired,float);f=np.asarray(faces,int)
    if baseline.shape!=target.shape or baseline.ndim!=3 or baseline.shape[-1]!=2:
        raise ValueError('Corresponding reference and target sections required')
    if not maximum_displacement_m>0 or not 0<minimum_relative_area<1:
        raise ValueError('Explicit positive displacement and relative-area limits required')
    area=lambda p:(p[f[:,1],0]-p[f[:,0],0])*(p[f[:,2],1]-p[f[:,0],1])-(p[f[:,1],1]-p[f[:,0],1])*(p[f[:,2],0]-p[f[:,0],0])
    result=target.copy();records=[];minimum=np.inf
    for k,(base,goal)in enumerate(zip(baseline,target)):
        original=area(base)
        if np.any(original<=0):raise ValueError('The reference section must already be strictly positive')
        delta=goal-base;lower,upper=0.,1.
        if np.any(area(goal)<original*minimum_relative_area):
            for _ in range(42):
                middle=(lower+upper)/2
                if np.all(area(base+delta*middle)>=original*minimum_relative_area):lower=middle
                else:upper=middle
            result[k]=base+delta*lower
        else:lower=1.
        distance=float(np.linalg.norm(result[k]-goal,axis=1).max())
        minimum=min(minimum,float((area(result[k])/original).min()))
        records.append(dict(station=k,retained_target_fraction=lower,maximum_target_deviation_m=distance))
    evidence=dict(passed=all(r['maximum_target_deviation_m']<=maximum_displacement_m for r in records),
        maximum_permitted_target_deviation_m=maximum_displacement_m,
        maximum_target_deviation_m=max(r['maximum_target_deviation_m']for r in records),
        minimum_relative_area=minimum,declared_minimum_relative_area=minimum_relative_area,
        adjusted_stations=[r for r in records if r['retained_target_fraction']<1.],
        source_envelope_and_anatomical_review_pending=True)
    return result,evidence


def transverse_segments(triangles,y):
    t=np.asarray(triangles,float);yy=t[:,:,1];t=t[(yy.min(1)<=y)&(yy.max(1)>=y)];result=[]
    for q in t:
        delta=q[:,1]-y;hits=[row[[0,2]]for row in q[delta==0]]
        for a,b in [(0,1),(1,2),(2,0)]:
            if delta[a]*delta[b]<0:hits.append((q[a]+(q[b]-q[a])*delta[a]/(delta[a]-delta[b]))[[0,2]])
        hits=np.unique(np.asarray(hits).reshape(-1,2),axis=0)
        if len(hits)==2:result.append(hits)
    if not result:raise ValueError('No transverse segments at the declared station')
    return np.asarray(result)


def measured_radial_curves(triangle_families,y,angles,*,duplicate_crossing_tolerance_m=2e-10):
    segments={k:transverse_segments(t,y)for k,t in triangle_families.items()};p=segments['inner'].reshape(-1,2)
    low,high=p.min(0),p.max(0);center=(low+high)/2;radii=(high-low)/2
    if np.any(radii<=0):raise ValueError('A section must have positive transverse dimensions')
    curves={};cross=lambda a,b:a[...,0]*b[...,1]-a[...,1]*b[...,0]
    for kind,s in segments.items():
        a=s[:,0]-center;b=s[:,1]-s[:,0];points=[]
        for theta in angles:
            direction=np.array([np.cos(theta),np.sin(theta)]);den=cross(direction,b);good=abs(den)>1e-20
            distance=np.divide(cross(a,b),den,out=np.full(len(a),np.nan),where=good)
            fraction=np.divide(cross(a,direction),den,out=np.full(len(a),np.nan),where=good)
            valid=good&(distance>0)&(fraction>=-1e-10)&(fraction<=1+1e-10)
            hits=distance[valid]
            if not len(hits):raise ValueError('No positive radial crossing')
            if hits.max()-hits.min()>duplicate_crossing_tolerance_m:raise ValueError('Multiple positive crossings: radial transport would change the occupied section')
            points.append(center+hits.mean()*direction)
        curves[kind]=np.asarray(points)
    if np.any(np.linalg.norm(curves['outer']-center,axis=1)<=np.linalg.norm(curves['inner']-center,axis=1)):
        raise ValueError('Pial inner/outer radial ordering failed')
    return center,radii,curves


def harmonic_sections(reference_positions,faces,fixed_controls,fixed_targets,*,transport_displacements=False,weights='uniform'):
    uv=np.asarray(reference_positions,float);f=np.asarray(faces,int);fixed=np.asarray(fixed_controls,int);targets=np.asarray(fixed_targets,float)
    if targets.ndim!=3 or targets.shape[1:]!=(len(fixed),2):raise ValueError('One corresponding two-dimensional fixed target per station required')
    if weights=='uniform':
        edges=np.unique(np.sort(np.r_[f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]],axis=1),axis=0)
        degree=np.bincount(edges.ravel(),minlength=len(uv));row=np.r_[np.arange(len(uv)),edges[:,0],edges[:,1]]
        col=np.r_[np.arange(len(uv)),edges[:,1],edges[:,0]];values=np.r_[degree,-np.ones(len(edges)*2)]
        laplacian=coo_matrix((values,(row,col)),shape=(len(uv),len(uv))).tocsc()
    elif weights=='mean-value':
        # Positive one-ring mean-value weights have affine precision, unlike
        # uniform graph weights on a nonuniform reference triangulation.
        # Floater 2003, https://doi.org/10.1016/S0167-8396(03)00002-5
        rows=[];cols=[];vals=[]
        for i,j,k in [(0,1,2),(1,2,0),(2,0,1)]:
            a=uv[f[:,j]]-uv[f[:,i]];b=uv[f[:,k]]-uv[f[:,i]]
            cross=abs(a[:,0]*b[:,1]-a[:,1]*b[:,0]);dot=np.einsum('ij,ij->i',a,b)
            half=np.tan(.5*np.arctan2(cross,dot))
            for vertex,edge in [(j,a),(k,b)]:
                rows.append(f[:,i]);cols.append(f[:,vertex]);vals.append(half/np.linalg.norm(edge,axis=1))
        matrix=coo_matrix((np.concatenate(vals),(np.concatenate(rows),np.concatenate(cols))),shape=(len(uv),len(uv))).tocsr()
        degree=np.asarray(matrix.sum(1)).ravel();assert np.all(degree>0)
        from scipy.sparse import diags
        laplacian=(diags(np.ones(len(uv)))-diags(1/degree)@matrix).tocsc()
    else:raise ValueError('Unknown explicitly requested reference weight scheme')
    interior=np.setdiff1d(np.arange(len(uv)),fixed)
    result=np.empty((len(targets),len(uv),2),float);result[:,fixed]=targets
    values=targets-uv[fixed] if transport_displacements else targets
    rhs=-(laplacian[interior][:,fixed]@values.transpose(1,0,2).reshape(len(fixed),-1))
    solved=splu(laplacian[interior][:,interior]).solve(rhs)
    result[:,interior]=solved.reshape(len(interior),len(targets),2).transpose(1,0,2)
    if transport_displacements:result[:,interior]+=uv[interior]
    q=result[:,f];area=(q[:,:,1,0]-q[:,:,0,0])*(q[:,:,2,1]-q[:,:,0,1])-(q[:,:,1,1]-q[:,:,0,1])*(q[:,:,2,0]-q[:,:,0,0])
    if np.any(area<=0):
        error=ValueError('Coupled harmonic transport inverted a reference section')
        error.bad_section_faces=np.argwhere(area<=0)
        error.section_fields=result
        error.minimum_section_double_area_m2=float(area.min())
        raise error
    return result,dict(stations=len(targets),controls_per_station=len(uv),minimum_double_area_m2=float(area.min()),
        all_triangles_positive=True,all_materials_share_one_map=True,topology_unchanged=True,
        transports_reference_displacements=transport_displacements,reference_weights=weights)


def transport_section_sequence(reference_positions,faces,fixed_controls,fixed_targets,anchor,*,interior_displacement_m=0.,wall_seconds=30.):
    """Transport nearby sections successively, rather than one distant warp.

    Positive mean-value weights are rebuilt from the last accepted geometry.
    An invalid step stops with its exact global station/triangle witnesses.
    """
    targets=np.asarray(fixed_targets,float);reference=np.asarray(reference_positions,float)
    if not 0<=anchor<len(targets):raise ValueError('An explicit in-range reference station is required')
    import time
    started=time.monotonic();result=np.full((len(targets),len(reference),2),np.nan);minimum=np.inf;repairs=[]
    for indices in [[anchor],list(range(anchor-1,-1,-1)),list(range(anchor+1,len(targets)))]:
        previous=reference if indices==[anchor]else result[anchor]
        for k in indices:
            if time.monotonic()-started>wall_seconds:raise TimeoutError('Successive section transport reached its native wall deadline')
            try:mapped,proof=harmonic_sections(previous,faces,fixed_controls,targets[k:k+1],weights='mean-value')
            except ValueError as error:
                if interior_displacement_m and hasattr(error,'section_fields'):
                    from geometry_tools.algorithms.section_kernel import repair_section_kernel,repair_section_patch
                    t=previous[np.asarray(faces,int)];a=t[:,1]-t[:,0];b=t[:,2]-t[:,0]
                    original_area=a[:,0]*b[:,1]-a[:,1]*b[:,0]
                    repaired,repair=repair_section_kernel(error.section_fields[0],faces,fixed_controls,
                        minimum_double_area=original_area*.05,maximum_displacement_m=interior_displacement_m)
                    repairs.append(dict(station=k,**repair))
                    if not repair['passed']:
                        repaired,repair=repair_section_patch(error.section_fields[0],faces,fixed_controls,
                            minimum_double_area=original_area*.05,maximum_displacement_m=interior_displacement_m,
                            initial=repaired)
                        repairs.append(dict(station=k,**repair))
                    if repair['passed']:
                        result[k]=repaired;previous=repaired;minimum=min(minimum,repair['minimum_double_area_m2'])
                        continue
                if hasattr(error,'bad_section_faces'):
                    result[k]=error.section_fields[0]
                    error.bad_section_faces[:,0]=k
                    error.section_fields=result
                    error.failed_sequence_station=k;error.sequence_repairs=repairs
                raise
            result[k]=mapped[0];minimum=min(minimum,proof['minimum_double_area_m2']);previous=result[k]
    return result,dict(stations=len(targets),controls_per_station=len(reference),minimum_double_area_m2=minimum,
        all_triangles_positive=True,all_materials_share_one_map=True,topology_unchanged=True,
        reference_weights='successive-mean-value',reference_station=anchor,interior_kernel_repairs=repairs,
        declared_interior_displacement_m=interior_displacement_m)


def radial_carrier_section(template, inner, outer, center, grey_faces, desired_grey_fraction,*,return_carrier=False):
    """Embed a conforming tissue graph through one positive affine carrier.

    Only the enclosing contours need be star shaped. The H-shaped grey,
    microscopic territories and central void are all carried by the same map.
    Radial grading changes the reference grey proportion without forcing its
    horns through the measured white-matter envelope. It is an anatomical
    estimate, not a claim of donor histological registration.
    """
    inner=np.asarray(inner,float);outer=np.asarray(outer,float);center=np.asarray(center,float)
    fractions=template['radial_fractions'];cf=template['carrier_faces'];owner=template['carrier_owner']
    bary=template['carrier_barycentric'];faces=template['faces'];n=len(inner)
    def cross(a,b):return a[...,0]*b[...,1]-a[...,1]*b[...,0]
    inner_area=.5*cross(inner-center,np.roll(inner,-1,axis=0)-center).sum()
    if inner_area<=0:raise ValueError('A measured inner contour must have positive orientation')
    def mapped(exponent):
        graded=np.where(fractions<=.1,fractions,.1+.9*(np.maximum(fractions-.1,0)/.9)**exponent)
        carrier=np.vstack([center[None],*[center+t*(inner-center)for t in graded],outer])
        q=carrier[cf];ca=cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
        if np.any(ca<=0):raise ValueError('The measured radial carrier has an inverted element')
        points=np.einsum('ni,nij->nj',bary,carrier[cf[owner]])
        q=points[faces];area=cross(q[:,1]-q[:,0],q[:,2]-q[:,0])
        return points,carrier,area,float(.5*area[grey_faces].sum()/inner_area),float(ca.min())
    low,high=.25,2.;a=mapped(low);b=mapped(high)
    if not b[3]<=desired_grey_fraction<=a[3]:
        raise ValueError('The explicit radial grading interval does not cover the requested regional grey proportion: '+
            str(dict(requested=desired_grey_fraction,range=[b[3],a[3]])))
    for _ in range(24):
        middle=(low+high)/2;trial=mapped(middle)
        if trial[3]>desired_grey_fraction:low=middle
        else:high=middle
    exponent=(low+high)/2;points,carrier,area,proportion,minimum=mapped(exponent)
    if np.any(area<=0):raise ValueError('A fine tissue triangle is not conforming to its positive carrier')
    # Fine triangles must also belong to a single carrier, not just have
    # positive vertices after unrelated point interpolation.
    tri_owner=template['triangle_carrier'];ref=template['carrier_positions'][cf[tri_owner]]
    mat=np.stack([ref[:,1]-ref[:,0],ref[:,2]-ref[:,0]],axis=2)
    weights=np.linalg.solve(mat[:,None],(template['positions'][faces]-ref[:,None,0])[...,None])[...,0]
    if np.min(np.c_[weights.reshape(-1,2),1-weights.sum(2).ravel()])< -1e-8:
        raise ValueError('A fine triangle crosses its declared affine carrier')
    proof=dict(all_triangles_positive=True,minimum_double_area_m2=float(area.min()),
        minimum_carrier_double_area_m2=minimum,radial_exponent=exponent,
        requested_grey_fraction=float(desired_grey_fraction),actual_grey_fraction=proportion,
        measured_inner_and_outer_controls_literal=True,central_void_is_transported=True,
        common_piecewise_affine_map=True)
    return (points,proof,carrier)if return_carrier else(points,proof)


def carrier_prism_edge_patches(barycentric, lower, upper):
    """Map one embedded edge through the three affine pieces of a prism.

    Barycentric coordinates use the globally sorted carrier vertices. The
    two internal planes are t=lambda_2 and t=1-lambda_0. Cutting at both gives
    exact shared longitudinal faces even for extremely narrow tissue regions;
    straight chords across different volume maps would not be equivalent.
    Returned polygons follow low-A, high-A, high-B, low-B winding.
    """
    cuts,points=carrier_prism_edge_grid(np.asarray(barycentric)[None],np.asarray(lower)[None],np.asarray(upper)[None])
    cuts=cuts[0];points=points[0]
    result=[]
    for piece in range(3):
        polygon=[]
        for endpoint,cut in [(0,piece),(0,piece+1),(1,piece+1),(1,piece)]:
            polygon.append((endpoint,float(cuts[endpoint,cut]),points[endpoint,cut]))
        result.append(polygon)
    return result


def carrier_prism_edge_grid(barycentric,lower,upper):
    """Batch the exact three-prism-map break points of embedded edges."""
    lam=np.asarray(barycentric,float);lower=np.asarray(lower,float);upper=np.asarray(upper,float)
    if lam.ndim!=3 or lam.shape[1:]!=(2,3)or lower.shape!=(len(lam),3,3)or upper.shape!=lower.shape:
        raise ValueError('Corresponding batches of edges and 3-D carrier prism triangles required')
    if lam.min()<-1e-8 or abs(lam.sum(2)-1).max()>1e-8:raise ValueError('The embedded edge leaves its declared carrier')
    lam=np.clip(lam,0,1);lam/=lam.sum(2)[...,None]
    a,b,c=[lam[:,:,i,None]for i in range(3)]
    cuts=np.stack([np.zeros(lam.shape[:2]),lam[:,:,2],1-lam[:,:,0],np.ones(lam.shape[:2])],axis=2)
    points=np.stack([np.einsum('eij,ejk->eik',lam,lower),
        a*lower[:,None,0]+b*lower[:,None,1]+c*upper[:,None,2],
        a*lower[:,None,0]+b*upper[:,None,1]+c*upper[:,None,2],
        np.einsum('eij,ejk->eik',lam,upper)],axis=2)
    return cuts,points


def prism_tetrahedra(fields,y,faces):
    fields=np.asarray(fields,float);y=np.asarray(y,float);f=np.sort(np.asarray(faces,int),axis=1)
    if len(fields)!=len(y)or np.any(np.diff(y)<=0):raise ValueError('Ordered distinct axial stations required')
    points=np.stack([np.c_[p[:,0],np.full(len(p),level),p[:,1]]for p,level in zip(fields,y)]).reshape(-1,3)
    n=fields.shape[1];tets=[];source=[]
    for k in range(len(y)-1):
        a,b,c=(f+k*n).T;d,e,h=(f+(k+1)*n).T
        tets.extend([np.c_[a,b,c,h],np.c_[a,b,e,h],np.c_[a,d,e,h]])
        source.extend([np.arange(len(f))]*3)
    tetrahedra=np.concatenate(tets);source=np.concatenate(source)
    det=np.linalg.det(points[tetrahedra[:,1:]]-points[tetrahedra[:,:1]])
    if np.any(det==0):raise ValueError('A calculation tetrahedron collapsed')
    first=fields[0,f];area=(first[:,1,0]-first[:,0,0])*(first[:,2,1]-first[:,0,1])-(first[:,1,1]-first[:,0,1])*(first[:,2,0]-first[:,0,0])
    expected=np.tile(np.r_[-np.sign(area),np.sign(area),-np.sign(area)],len(y)-1)
    wrong=np.flatnonzero(np.sign(det)!=expected)
    if len(wrong):raise PrismJacobianError(np.unique(wrong//(3*len(f))))
    negative=expected<0;tetrahedra[negative]=tetrahedra[negative][:,[0,2,1,3]]
    # Ordering the face controls makes every shared prism quadrilateral choose
    # the same diagonal. Signed volume is independently retained for replay.
    return points,tetrahedra,source,dict(tetrahedra=len(tetrahedra),all_positive_reference_jacobians=True,
        minimum_tetrahedron_volume_m3=float(abs(det).min()/6),volume_m3=float(abs(det).sum()/6))
