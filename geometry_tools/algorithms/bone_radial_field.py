"""Measured radial placement ratios for corresponding fixed digit bones.

This moves section centres only. It does not resize a vessel, change a bone,
certify containment, or infer anatomical ownership. Ray misses stay recorded;
nearest valid cylindrical samples only initialize an explicitly provisional fit.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator,RegularGridInterpolator
from scipy.ndimage import distance_transform_edt,gaussian_filter
from geometry_tools.open3d_backend import SurfaceIndex,library
from geometry_tools.oriented_rays import literal_hits
from .tube_sections import rotation_minimizing_frames,unit


def last_surface_hit_before(index,origins,directions,limits_m):
    """Last literal surface hit on each finite forward ray, before its limit.

    Zero means no validated hit on this finite ray, not empty solid volume.
    Open sheets, tissue identity and complete-profile clearance need separate
    evidence. Returned face IDs permit independent material-owner replay.
    """
    origins=np.asarray(origins,float);directions=np.asarray(directions,float);limits=np.asarray(limits_m,float)
    if origins.ndim!=2 or origins.shape[1]!=3 or directions.shape!=origins.shape or limits.shape!=(len(origins),):raise ValueError('One origin, unit direction and positive finite limit per ray required')
    if not np.isfinite(np.r_[origins.ravel(),directions.ravel(),limits]).all()or np.any(limits<=0)or not np.allclose(np.linalg.norm(directions,axis=1),1,atol=1e-10,rtol=0):raise ValueError('Finite forward ray intervals and unit directions required')
    rays=np.c_[(origins-index.origin)/index.scale,directions].astype(np.float32)
    raw=index.scene.list_intersections(library().core.Tensor(rays),nthreads=2)
    ids=raw['ray_ids'].numpy().astype(int);faces=raw['primitive_ids'].numpy().astype(int)
    valid,travel,_,_,_=literal_hits(origins[ids],directions[ids],index.mesh.triangles[faces])
    valid&=(travel>1e-9)&(travel<limits[ids]);radii=np.zeros(len(origins));np.maximum.at(radii,ids[valid],travel[valid])
    selected=valid&(travel==radii[ids]);owner=np.full(len(origins),len(index.mesh.faces),int);np.minimum.at(owner,ids[selected],faces[selected]);owner[owner==len(index.mesh.faces)]=-1
    return dict(distance_m=radii,face_ids=owner,validated_hits=int(valid.sum()),no_hit_rays=int((owner<0).sum()),solidClearanceCertified=False)


def _outer_radii(index,origins,directions,maximum_radius_m,*,first=False):
    rays=np.c_[(origins-index.origin)/index.scale,directions].astype(np.float32)
    raw=index.scene.list_intersections(library().core.Tensor(rays),nthreads=2)
    ids=raw['ray_ids'].numpy().astype(int);faces=raw['primitive_ids'].numpy().astype(int)
    valid,travel,_,_,_=literal_hits(origins[ids],directions[ids],index.mesh.triangles[faces])
    valid&=(travel<maximum_radius_m)
    radii=np.full(len(origins),np.inf)if first else np.zeros(len(origins))
    (np.minimum.at if first else np.maximum.at)(radii,ids[valid],travel[valid])
    radii[~np.isfinite(radii)]=0.
    return radii


def cached_bone_radial_field(bank):
    """Replay a saved measurement grid; caller must bind its source/hash proof."""
    t=np.asarray(bank['stations'],float);angles=np.asarray(bank['angles'],float)
    frames=np.asarray(bank['source_frames'],float);radii=np.asarray(bank['annulus_radii'],float)
    if t.ndim!=1 or np.any(np.diff(t)<=0)or angles.ndim!=1 or np.any(np.diff(angles)<=0):
        raise ValueError('Ordered saved stations and periodic angles required')
    if frames.shape!=(len(t),3,3)or radii.shape!=(4,len(t),len(angles)):
        raise ValueError('Corresponding saved frame and bone/skin measurement grids required')
    if not np.isfinite(radii).all()or np.any(radii<=0)or np.any(radii[1]<=radii[0])or np.any(radii[3]<=radii[2]):
        raise ValueError('Positive ordered bone/skin intervals required')
    if not np.allclose(frames.transpose(0,2,1)@frames,np.eye(3),atol=1e-8,rtol=0)or not np.allclose(np.linalg.det(frames),1,atol=1e-8,rtol=0):
        raise ValueError('Saved source frames must be proper rotations')
    if abs(angles[0])>1e-12 or abs(angles[-1]-2*np.pi)>1e-12 or not np.array_equal(radii[:,:,-1],radii[:,:,0]):
        raise ValueError('Literal periodic angular seam required')
    interp=[RegularGridInterpolator((t,angles),r,bounds_error=True)for r in radii]
    def field(station,delta):
        station=np.clip(np.asarray(station,float),t[0],t[-1]);delta=np.asarray(delta,float)
        k=np.clip(np.searchsorted(t,station)-1,0,len(t)-2);u=(station-t[k])/(t[k+1]-t[k])
        frame=(1-u[:,None,None])*frames[k]+u[:,None,None]*frames[k+1]
        angle=np.mod(np.arctan2(np.einsum('ni,ni->n',delta,frame[:,:,1]),np.einsum('ni,ni->n',delta,frame[:,:,0])),2*np.pi)
        sb,ss,tb,ts=[f(np.c_[station,angle])for f in interp];radius=np.linalg.norm(delta,axis=1)
        mapped=tb+(radius-sb)*(ts-tb)/(ss-sb);inside=radius<sb
        mapped[inside]=radius[inside]*tb[inside]/sb[inside]
        return np.divide(mapped,radius,out=tb/sb,where=radius>1e-12)
    return field


def measured_bone_radial_field(source_mesh,target_mesh,source_bones,target_bones,
                               source_normal,target_normal,*,maximum_radius_m=.025,
                               source_skin=None,target_skin=None):
    a=np.asarray(source_bones,float);b=np.asarray(target_bones,float)
    if a.shape!=b.shape or a.shape[1:]!=(3,)or len(a)<3:
        raise ValueError('Explicit corresponding bone centres required')
    if not .01<=maximum_radius_m<=.03:raise ValueError('Explicit local finger radius bound required')
    t=np.linspace(-.75,len(a)-.25,257);angles=np.linspace(0,2*np.pi,96,endpoint=False)
    if (source_skin is None)!=(target_skin is None):raise ValueError('Both source and target skin envelopes required')
    frames=[];radii=[];skin_radii=[]
    for mesh,bones,normal,skin in [(source_mesh,a,source_normal,source_skin),(target_mesh,b,target_normal,target_skin)]:
        controls=np.r_[(2*bones[0]-bones[1])[None],bones,(2*bones[-1]-bones[-2])[None]]
        curve=PchipInterpolator(np.arange(-1,len(bones)+1),controls);axis=curve(t)
        frame=rotation_minimizing_frames(axis,unit(curve.derivative()(t)),initial_reference=normal)
        directions=frame[:,None,:,0]*np.cos(angles)[None,:,None]+frame[:,None,:,1]*np.sin(angles)[None,:,None]
        origins=np.broadcast_to(axis[:,None],directions.shape)
        r=_outer_radii(SurfaceIndex(mesh),origins.reshape(-1,3),directions.reshape(-1,3),maximum_radius_m).reshape(len(t),len(angles))
        frames.append(frame);radii.append(r)
        if skin is not None:
            skin_radii.append(_outer_radii(SurfaceIndex(skin),origins.reshape(-1,3),directions.reshape(-1,3),.05,first=True).reshape(len(t),len(angles)))
    valid=(radii[0]>.001)&(radii[1]>.001)
    if not valid.any():raise ValueError('No corresponding bone ray samples')
    ratio=np.divide(radii[1],radii[0],out=np.full_like(radii[0],np.nan),where=valid)
    # The finite ratio bounds restrict this initializer; clipped samples are
    # reported as unresolved rather than silently accepted as real anatomy.
    clipped=valid&((ratio<.5)|(ratio>3.))
    ratio[valid]=np.clip(ratio[valid],.5,3.)
    # Periodic angular padding prevents a false discontinuity at atan2=0.
    padded=np.tile(ratio,(1,3));known=np.isfinite(padded)
    nearest=distance_transform_edt(~known,return_distances=False,return_indices=True)
    filled=padded[tuple(nearest)]
    smooth=gaussian_filter(filled,sigma=(2.,1.5),mode=('nearest','wrap'))[:,96:192]
    interp=RegularGridInterpolator((t,np.r_[angles,2*np.pi]),np.c_[smooth,smooth[:,0]],bounds_error=True)
    annulus=[];annulus_samples=[];interval_evidence=None
    if skin_radii:
        interval_valid=valid&(skin_radii[0]>radii[0]+.0005)&(skin_radii[1]>radii[1]+.0005)
        if not interval_valid.any():raise ValueError('No corresponding finite bone/skin corridors')
        for values in [radii[0],skin_radii[0],radii[1],skin_radii[1]]:
            pad=np.tile(np.where(interval_valid,values,np.nan),(1,3));known=np.isfinite(pad)
            nearest=distance_transform_edt(~known,return_distances=False,return_indices=True)
            smooth_values=gaussian_filter(pad[tuple(nearest)],sigma=(2.,1.5),mode=('nearest','wrap'))[:,96:192]
            values=np.c_[smooth_values,smooth_values[:,0]];annulus_samples.append(values)
            annulus.append(RegularGridInterpolator((t,np.r_[angles,2*np.pi]),values,bounds_error=True))
        interval_evidence=dict(measuredPairedCorridors=int(interval_valid.sum()),unmeasuredPairedCorridors=int((~interval_valid).sum()),
            sourceSkinSha256=source_skin.geometry_sha,targetSkinSha256=target_skin.geometry_sha,
            sourceMinimumMeasuredWidth_m=float((skin_radii[0]-radii[0])[interval_valid].min()),
            targetMinimumMeasuredWidth_m=float((skin_radii[1]-radii[1])[interval_valid].min()),
            mapping='Piecewise monotone radius: bone interior scales from the axis; soft tissue preserves normalized bone/skin position; section gauge retained',
            interiorBoneQuerySamples=0,exteriorSkinQuerySamples=0,
            queriesIncludeRawControlsAndSectionCentres=True,sourceBoneEntryAnatomyNotInferred=True)
    def field(station,delta):
        station=np.clip(station,t[0],t[-1]);k=np.clip(np.searchsorted(t,station)-1,0,len(t)-2)
        u=(station-t[k])/(t[k+1]-t[k])
        frame=(1-u[:,None,None])*frames[0][k]+u[:,None,None]*frames[0][k+1]
        angle=np.mod(np.arctan2(np.einsum('ni,ni->n',delta,frame[:,:,1]),np.einsum('ni,ni->n',delta,frame[:,:,0])),2*np.pi)
        coords=np.c_[station,angle]
        if annulus:
            sb,ss,tb,ts=[f(coords)for f in annulus];radius=np.linalg.norm(delta,axis=1)
            fraction=(radius-sb)/(ss-sb)
            mapped=tb+fraction*(ts-tb)
            interior=radius<sb
            interval_evidence['interiorBoneQuerySamples']+=int(interior.sum())
            interval_evidence['exteriorSkinQuerySamples']+=int((radius>ss).sum())
            # Extending the annular affine map through the axis can reverse
            # radius. A positive inner scaling instead preserves the donor's
            # intracortical relation without deleting or relabeling a branch.
            mapped[interior]=radius[interior]*tb[interior]/sb[interior]
            scale=np.divide(mapped,radius,out=tb/sb,where=radius>1e-12)
            if np.any(scale<=0):raise ValueError('Nonpositive measured radial continuation')
            return scale
        return interp(coords)
    proof=dict(method='Literal Float64 replay of native all-hit proposals on corresponding fixed bone envelopes',
        samples=int(valid.size),correspondingMeasuredSamples=int(valid.sum()),unmeasuredSamples=int((~valid).sum()),
        clippedInitializerSamples=int(clipped.sum()),measuredRatioMedian=float(np.nanmedian(ratio)),
        proposedRatioMinimum=float(smooth.min()),proposedRatioMaximum=float(smooth.max()),
        sourceGeometrySha256=source_mesh.geometry_sha,targetGeometrySha256=target_mesh.geometry_sha,
        maximumRayRadius_m=maximum_radius_m,unmeasuredSamplesInterpolatedForProposalOnly=True,
        sourceSectionCalibreChanged=False,skeletonChanged=False,anatomicalAcceptance=False)
    proof['pairedBoneSkinIntervals']=interval_evidence
    field.sample_bank=dict(stations=t,angles=np.r_[angles,2*np.pi],source_frames=frames[0],
        radial_ratios=np.c_[smooth,smooth[:,0]],annulus_radii=np.array(annulus_samples))
    return field,proof
