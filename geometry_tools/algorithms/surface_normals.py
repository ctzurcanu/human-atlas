"""Smooth explicitly grouped surface charts across exact coordinate aliases.

The caller supplies a group for every triangle. No tissue identity, crease
angle, proximity weld or anatomical continuity is inferred by this function.
"""
import numpy as np


def transport_vertex_words(reference_positions, positions, faces, words, *, maximum_facet_rotation_degrees):
    """Transport a custom field through corresponding affine surface charts.

    Each triangle supplies two actual edge derivatives and a unit normal
    extension. Source-area weighted inverse-transpose corner fields are shared
    only at exact coordinate aliases. Unchanged stars retain literal words.
    This changes shading only, never positions, source topology or UVs.
    """
    old=np.asarray(reference_positions,float);new=np.asarray(positions,float)
    f=np.asarray(faces,int);w=np.asarray(words)
    if old.shape!=new.shape or old.shape!=w.shape or not np.isfinite(np.r_[old.ravel(),new.ravel(),w.ravel()]).all():
        raise ValueError('Finite corresponding source controls and normal words required')
    bound=float(maximum_facet_rotation_degrees)
    if not 0<=bound<90:raise ValueError('Explicit facet rotation bound below 90 degrees required')
    a=old[f];b=new[f];ea=a[:,1:]-a[:,:1];eb=b[:,1:]-b[:,:1]
    na=np.cross(ea[:,0],ea[:,1]);nb=np.cross(eb[:,0],eb[:,1])
    la=np.linalg.norm(na,axis=1);lb=np.linalg.norm(nb,axis=1)
    if np.any(la==0)or np.any(lb==0):raise ValueError('Nondegenerate corresponding facets required')
    na/=la[:,None];nb/=lb[:,None]
    angles=np.degrees(np.arccos(np.clip(np.einsum('ij,ij->i',na,nb),-1,1)))
    if angles.max(initial=0)>bound+1e-9:raise ValueError('Facet rotation exceeds declared transport bound')
    A=np.stack([ea[:,0],ea[:,1],na],axis=2);B=np.stack([eb[:,0],eb[:,1],nb],axis=2)
    J=np.linalg.solve(A.transpose(0,2,1),B.transpose(0,2,1)).transpose(0,2,1)
    field=np.linalg.solve(J.transpose(0,2,1),w[f].astype(float).transpose(0,2,1)).transpose(0,2,1)
    lengths=np.linalg.norm(field,axis=2)
    if np.any(lengths==0):raise ValueError('Source field contains zero normals')
    field/=lengths[:,:,None]
    _,inverse=np.unique(old,axis=0,return_inverse=True);n=int(inverse.max())+1
    sums=np.zeros((n,3));active=np.zeros(n,bool)
    displacement=new-old;changed=np.any(displacement[f]!=displacement[f[:,:1]],axis=(1,2))
    for k in range(3):
        np.add.at(sums,inverse[f[:,k]],field[:,k]*la[:,None])
        np.logical_or.at(active,inverse[f[:,k]],changed)
    length=np.linalg.norm(sums,axis=1)
    if np.any(active&(length==0)):raise ValueError('Transported incident fields cancel')
    result=w.astype('<i2',copy=True);take=active[inverse]
    result[take]=np.rint(sums[inverse[take]]/length[inverse[take],None]*32767).astype('<i2')
    return result,dict(method='source-area-weighted affine chart inverse-transpose',
        maximum_facet_rotation_degrees=float(angles.max(initial=0)),changed_vertex_fields=int(take.sum()),
        unchanged_stars_literal=bool(np.array_equal(result[~take],w[~take])),positions_modified=False)


def transport_corner_words(reference_positions,positions,faces,words,*,maximum_facet_rotation_degrees):
    """Carry a corner field by each facet's proper normal rotation.

    Use only for an explicitly deforming chart. This preserves each corner's
    angle to its facet before the small outward-cone encoding correction;
    reciprocal material owners still reuse the same words with opposite sign.
    """
    reference=np.asarray(reference_positions,float);p=np.asarray(positions,float);f=np.asarray(faces,int);w=np.asarray(words,float)
    bound=float(maximum_facet_rotation_degrees)
    if reference.shape!=p.shape or w.shape!=(len(f),3,3)or not 0<=bound<90:raise ValueError('Corresponding charts and an explicit rotation bound below 90 degrees required')
    old=reference[f];new=p[f];a=np.cross(old[:,1]-old[:,0],old[:,2]-old[:,0]);b=np.cross(new[:,1]-new[:,0],new[:,2]-new[:,0]);la=np.linalg.norm(a,axis=1);lb=np.linalg.norm(b,axis=1)
    if np.any(la==0)or np.any(lb==0)or not np.isfinite(np.r_[a.ravel(),b.ravel(),w.ravel()]).all():raise ValueError('Finite nondegenerate charts required')
    a/=la[:,None];b/=lb[:,None];cos=np.clip(np.einsum('fi,fi->f',a,b),-1,1);angle=np.degrees(np.arccos(cos))
    if angle.max(initial=0)>bound+1e-9:raise ValueError('Facet rotation exceeds the declared normal-transport bound')
    if np.any(np.einsum('fki,fi->fk',w,a)<=0):raise ValueError('Source corner field must be outward')
    axis=np.cross(a,b)[:,None];unit=w/np.linalg.norm(w,axis=2)[:,:,None]
    rotated=unit+np.cross(axis,unit)+np.cross(axis,np.cross(axis,unit))/(1+cos[:,None,None])
    encoded,evidence=quantize_outward_corner_words(p,f,rotated,maximum_adjustment_degrees=.3,minimum_outward_cosine=.002)
    actual=encoded.astype(float);actual/=np.linalg.norm(actual,axis=2)[:,:,None]
    directions=np.degrees(np.arccos(np.clip(np.einsum('fki,fki->fk',actual,unit),-1,1)))
    return encoded,dict(facets=len(f),maximum_facet_rotation_degrees=float(angle.max(initial=0)),maximum_corner_rotation_degrees=float(directions.max(initial=0)),outward_encoding=evidence,positions_faces_uvs_unchanged=True)


def quantize_outward_corner_words(positions, faces, normals, *, maximum_adjustment_degrees,
                                 minimum_outward_cosine=0.002):
    """Encode a transported corner field with an explicit outward-angle bound.

    A positive anisotropic map need not preserve the acute angle between two
    independently transported normals. Only offending shading corners are
    projected to the requested outward cone. Positions, topology and UVs are
    untouched; this is neither a geometric repair nor an anatomical assertion.
    """
    p=np.asarray(positions,float); f=np.asarray(faces,int); n=np.asarray(normals,float)
    bound=float(maximum_adjustment_degrees); margin=float(minimum_outward_cosine)
    if n.shape!=(len(f),3,3) or not 0<margin<1 or not 0<=bound<90:
        raise ValueError('Corresponding corner normals and explicit finite angular bounds are required')
    if not np.isfinite(np.r_[p.ravel(),n.ravel(),bound,margin]).all():
        raise ValueError('Finite normal encoding inputs are required')
    q=p[f]; area=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]); length=np.linalg.norm(area,axis=1)
    source_length=np.linalg.norm(n,axis=2)
    if np.any(length==0) or np.any(source_length==0):
        raise ValueError('Nondegenerate facets and nonzero corner fields are required')
    facet=area/length[:,None]; source=n/source_length[:,:,None]; adjusted=source.copy()
    cosine=np.einsum('tki,ti->tk',source,facet)
    bad=cosine<margin
    tangential=source-cosine[:,:,None]*facet[:,None]
    tangent_length=np.linalg.norm(tangential,axis=2)
    if np.any(bad & (tangent_length==0)):
        raise ValueError('An antiparallel field exceeds the bounded outward correction')
    tangential=np.divide(tangential,tangent_length[:,:,None],out=np.zeros_like(tangential),
                         where=tangent_length[:,:,None]>0)
    projected=tangential*np.sqrt(1-margin**2)+facet[:,None]*margin
    adjusted[bad]=projected[bad]
    words=np.rint(adjusted*32767).astype('<i2')
    encoded=words.astype(float); encoded/=np.linalg.norm(encoded,axis=2)[:,:,None]
    angle=np.degrees(np.arccos(np.clip(np.einsum('tki,tki->tk',encoded,source),-1,1)))
    largest=float(angle.max(initial=0))
    if largest>bound or np.any(np.einsum('tki,ti->tk',encoded,facet)<=0):
        raise ValueError('The requested normal adjustment bound cannot encode an outward field')
    return words,dict(adjusted_corners=int(bad.sum()),maximum_adjustment_degrees=largest,
        declared_adjustment_bound_degrees=bound,minimum_outward_cosine=margin,
        all_corner_words_outward=True,positions_modified=False,topology_modified=False)


def transport_z_shear_words(positions, faces, vertex_words, displacement_z):
    """Transport a shading field by an explicit piecewise Z-shear gradient.

    Source-area-weighted incident XY gradients define the field at physical
    vertices. Exact coordinate aliases receive the same gradient. Unchanged
    stars retain their literal words. A changed vertical chart needs a different
    deformation model and is rejected instead of guessing its derivative.
    """
    p=np.asarray(positions,float);f=np.asarray(faces,int);words=np.asarray(vertex_words)
    d=np.asarray(displacement_z,float)
    if p.ndim!=2 or p.shape[1]!=3 or words.shape!=p.shape or d.shape!=(len(p),) \
            or not np.isfinite(np.r_[p.ravel(),words.ravel(),d]).all():
        raise ValueError('Finite corresponding positions, normal words and Z displacements are required')
    t=p[f];delta=d[f];changed=np.any(delta!=delta[:,:1],axis=1)
    gradient=np.zeros((len(f),2))
    if np.any(changed):
        edges=t[changed,1:,:2]-t[changed,:1,:2]
        determinant=np.linalg.det(edges)
        if np.any(determinant==0):
            raise ValueError('A changed vertical chart has no affine XY shear gradient')
        rhs=delta[changed,1:]-delta[changed,:1]
        gradient[changed]=np.linalg.solve(edges,rhs[:,:,None])[:,:,0]
    area=np.linalg.norm(np.cross(t[:,1]-t[:,0],t[:,2]-t[:,0]),axis=1)/2
    if np.any(area==0):raise ValueError('A degenerate source face cannot transport a normal field')
    unique,inverse=np.unique(p,axis=0,return_inverse=True)
    total=np.zeros(len(unique));sums=np.zeros((len(unique),2))
    for k in range(3):
        np.add.at(total,inverse[f[:,k]],area)
        np.add.at(sums,inverse[f[:,k]],gradient*area[:,None])
    average=np.divide(sums,total[:,None],out=np.zeros_like(sums),where=total[:,None]>0)[inverse]
    active=np.any(average!=0,axis=1);result=words.astype('<i2',copy=True)
    n=words[active].astype(float);n[:,:2]-=average[active]*n[:,2,None]
    lengths=np.linalg.norm(n,axis=1)
    if np.any(lengths==0):raise ValueError('The supplied shading field contains a zero normal')
    result[active]=np.rint(np.clip(n/lengths[:,None],-1,1)*32767).astype('<i2')
    return result,dict(method='source-area-weighted affine XY shear gradients; inverse-transpose normal transport',
        changed_face_gradients=int(changed.sum()),changed_vertex_fields=int(active.sum()),
        unchanged_vertex_words_literal=bool(np.array_equal(result[~active],words[~active])),
        positions_modified=False)


def source_normal_words(positions,surfaces,parent_faces,source_positions,source_faces,source_corner_words):
    """Transfer one physical source field independently of tissue tessellation."""
    from geometry_tools.distance import barycentric
    v=np.asarray(positions,float);sv=np.asarray(source_positions,float);sf=np.asarray(source_faces,int)
    sn=np.asarray(source_corner_words,float);result={};fallbacks=0
    if sn.shape!=(len(sf),3,3)or set(surfaces)!=set(parent_faces):
        raise ValueError('Explicit source corner fields and parent facets are required')
    for sid,faces in surfaces.items():
        faces=np.asarray(faces,int);parents=np.asarray(parent_faces[sid],int)
        if parents.shape!=(len(faces),)or np.any((parents<0)|(parents>=len(sf))):
            raise ValueError('Every actual physical triangle requires a valid source facet')
        q=v[faces];triangles=sv[sf[parents]]
        weights,valid=barycentric(q.reshape(-1,3),np.repeat(triangles,3,axis=0))
        if not np.all(valid):raise ValueError('A source shading facet has zero area')
        weights=weights.reshape(-1,3,3)
        # Quantization can put a corner just outside its source face. Only the
        # shading sample is clamped; no geometric position is projected/moved.
        weights=np.clip(weights,0,1);weights/=weights.sum(2)[:,:,None]
        normals=np.einsum('tki,tij->tkj',weights,sn[parents])
        normals/=np.linalg.norm(normals,axis=2)[:,:,None]
        fn=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);length=np.linalg.norm(fn,axis=1)
        if np.any(length==0):raise ValueError('Stored physical triangles must have positive area')
        bad=np.einsum('tki,ti->tk',normals,fn)<=0
        normals[bad]=np.broadcast_to((fn/length[:,None])[:,None],normals.shape)[bad]
        words=np.rint(np.clip(normals,-1,1)*32767).astype('<i2')
        if np.any(np.einsum('tki,ti->tk',words.astype(float),fn)<=0):
            raise ValueError('Transferred physical source normals are not outward')
        result[sid]=words;fallbacks+=int(bad.sum())
    return result,dict(method='barycentric physical source corner normals; independent of tissue partitions',
        hard_corner_fallbacks=fallbacks,all_corner_words_outward=True,source_facet_witnesses=True,geometric_positions_modified=False)


def grouped_normal_words(positions, surfaces, face_groups):
    positions=np.asarray(positions,float)
    if set(surfaces)!=set(face_groups):
        raise ValueError('Every surface requires explicit triangle normal groups')
    vertices,inverse=np.unique(positions,axis=0,return_inverse=True)
    groups={};result={};fallbacks=0;quantization_fallbacks=0
    for sid,faces in surfaces.items():
        faces=np.asarray(faces,int);labels=np.asarray(face_groups[sid],object)
        if labels.shape!=(len(faces),):
            raise ValueError('A normal group must cover each actual surface triangle')
        result[sid]=np.empty((len(faces),3,3),dtype='<i2')
        for group in set(labels):
            ids=np.flatnonzero(labels==group)
            groups.setdefault(group,[]).append((sid,ids,inverse[faces[ids]]))
    for entries in groups.values():
        fs=np.concatenate([f for _,_,f in entries]);q=vertices[fs]
        normals=np.cross(q[:,1]-q[:,0],q[:,2]-q[:,0]);length=np.linalg.norm(normals,axis=1)
        if np.any(length==0):
            raise ValueError('A zero-area triangle cannot receive an outward normal')
        sums=np.zeros_like(vertices)
        for k in range(3):
            np.add.at(sums,fs[:,k],normals)
        corners=sums[fs];norm=np.linalg.norm(corners,axis=2)
        corners=np.divide(corners,norm[:,:,None],out=np.zeros_like(corners),where=norm[:,:,None]>0)
        bad=np.einsum('tki,ti->tk',corners,normals)<=0
        corners[bad]=np.broadcast_to((normals/length[:,None])[:,None],corners.shape)[bad]
        fallbacks+=int(bad.sum())
        words=np.rint(np.clip(corners,-1,1)*32767).astype('<i2')
        # A nearly tangent but outward Float64 average can round to a tangent
        # or inward Int16 vector. Retain the actual facet's outward direction
        # for only those corners, just as for the pre-encoding hard fallback.
        quantized_bad=np.einsum('tki,ti->tk',words.astype(float),normals)<=0
        if np.any(quantized_bad):
            facet=np.rint(np.clip(normals/length[:,None],-1,1)*32767).astype('<i2')
            words[quantized_bad]=np.broadcast_to(facet[:,None],words.shape)[quantized_bad]
            quantization_fallbacks+=int(quantized_bad.sum())
        if np.any(np.einsum('tki,ti->tk',words.astype(float),normals)<=0):
            raise ValueError('Quantized chart normals are not outward')
        offset=0
        for sid,ids,faces in entries:
            result[sid][ids]=words[offset:offset+len(faces)];offset+=len(faces)
    return result,dict(method='area-weighted explicit groups, exact coordinate aliases only',
                       groups=len(groups),original_controls=len(positions),unique_controls=len(vertices),
                       hard_corner_fallbacks=fallbacks,quantized_hard_corner_fallbacks=quantization_fallbacks,all_corner_words_outward=True)
