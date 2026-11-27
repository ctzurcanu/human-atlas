"""Transport source sections by corresponding named bone stations.

The paired skeletons are guides, never movable operands. A proper pose at
each original section preserves its actual calibre and branch ownership.
This does not establish a free vascular corridor or anatomical acceptance.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from .tube_axis_map import axis_transport
from .tube_sections import rotation_minimizing_frames,unit


def paired_curve_carriage(vertices,sections,physical_ids,source_curves,target_curves,
                          source_axis,target_axis,source_origin,target_origin,axial_scale):
    """Carry an intact network between two explicit corresponding curved rails.

    Rails define a shared transverse frame at each axial station; they are
    proposals, not certified free corridors. Centres follow the varying rail
    width. Every actual source section receives only a proper rigid pose.
    """
    from scipy.spatial.transform import Rotation,Slerp
    v=np.asarray(vertices,float);physical=np.asarray(physical_ids,int)
    if physical.shape!=(len(v),):raise ValueError('Original physical identities required')
    if not .5<=axial_scale<=2.5:raise ValueError('Finite positive bounded axial scale required')
    axes=[unit(source_axis),unit(target_axis)];origins=[np.asarray(source_origin,float),np.asarray(target_origin,float)]
    curves=[[np.asarray(c,float)for c in pair]for pair in [source_curves,target_curves]]
    if any(c.ndim!=2 or c.shape[1]!=3 or len(c)<3 or not np.isfinite(c).all() for pair in curves for c in pair):
        raise ValueError('Two finite ordered source and recipient courses required')
    parameters=[[((c-o)@axis)for c in pair]for pair,o,axis in zip(curves,origins,axes)]
    lo=max(min(s)for s in parameters[0]);hi=min(max(s)for s in parameters[0])
    if hi-lo<.001:raise ValueError('Source rail overlap too short')
    station=np.linspace(lo,hi,257)
    def sample(pair,params,values,axis,origin):
        result=[]
        for c,s in zip(pair,params):
            if np.min(np.diff(s))<-.001:raise ValueError('A rail reverses by more than 1 mm; explicit longitudinal correspondence required')
            order=np.argsort(s);s=s[order];c=c[order]
            if np.any(np.diff(s)<1e-12):raise ValueError('Repeated rail station requires explicit correspondence')
            q=np.column_stack([np.interp(values,s,c[:,j])for j in range(3)])
            # Preserve the requested plane even outside the measured rail range;
            # endpoint transverse offsets stay constant, never invent a bend.
            q+=(values-(q-origin)@axis)[:,None]*axis
            result.append(q)
        return result
    a=sample(curves[0],parameters[0],station,axes[0],origins[0])
    b=sample(curves[1],parameters[1],station*axial_scale,axes[1],origins[1])
    def frame(pair,axis):
        delta=pair[0]-pair[1];delta-=np.outer(delta@axis,axis);width=np.linalg.norm(delta,axis=1)
        x=delta/np.maximum(width[:,None],1e-30);z=np.broadcast_to(axis,x.shape)
        return np.stack([x,np.cross(z,x),z],axis=2),width
    before,wa=frame(a,axes[0]);after,wb=frame(b,axes[1])
    valid=(wa>.001)&(wb>.001)&(wb/np.maximum(wa,1e-30)>=.5)&(wb/np.maximum(wa,1e-30)<=2.5)
    ids=np.flatnonzero(valid)
    if len(ids)<3 or np.any(np.diff(ids)!=1):raise ValueError('No single separated bounded-width paired-rail interval')
    station=station[ids];before=before[ids];after=after[ids];width=wb[ids]/wa[ids]
    old=(a[0][ids]+a[1][ids])*.5;new=(b[0][ids]+b[1][ids])*.5
    def place(points):
        s=(points-origins[0])@axes[0];u=np.clip(s,station[0],station[-1])
        F=Slerp(station,Rotation.from_matrix(before))(u).as_matrix()
        G=Slerp(station,Rotation.from_matrix(after))(u).as_matrix()
        a0=np.column_stack([np.interp(u,station,old[:,j])for j in range(3)])
        b0=np.column_stack([np.interp(u,station,new[:,j])for j in range(3)])
        local=np.einsum('nji,nj->ni',F,points-a0)
        local[:,:2]*=np.interp(u,station,width)[:,None];local[:,2]*=axial_scale
        return b0+np.einsum('nij,nj->ni',G,local),G@F.transpose(0,2,1)
    q,R=place(v);owned=np.zeros(len(v),bool)
    for section in sections:
        s=np.asarray(section,int)
        if owned[s].any():raise ValueError('Original source sections overlap')
        centre=v[s].mean(0);point,rotation=place(centre[None])
        for raw in s:
            aliases=physical==physical[raw];q[aliases]=(v[aliases]-centre)@rotation[0].T+point[0];R[aliases]=rotation[0];owned[aliases]=True
    return q,R,dict(method='Two corresponding curved rails; shared axial planes and proper section poses',
                    sourceSectionCalibreUnscaled=True,pairedStationInterval_m=station[[0,-1]].tolist(),
                    transverseCarriageScaleRange=[float(width.min()),float(width.max())],axialScale=float(axial_scale),
                    railClearanceCertified=False,anatomicalAcceptance=False)


def paired_vascular_frame_carriage(vertices,sections,physical_ids,source_pair,target_pair,source_axis,target_axis,axial_scale):
    """A positive two-rail carrier; vessel cross-sections retain literal shape."""
    v=np.asarray(vertices,float);a,b=np.asarray(source_pair,float),np.asarray(target_pair,float)
    if a.shape!=(2,3)or b.shape!=(2,3)or not np.isfinite(np.r_[a.ravel(),b.ravel()]).all():
        raise ValueError('Two explicit corresponding vascular reference points per frame required')
    def frame(pair,axis):
        z=unit(axis);x=pair[0]-pair[1];x-=z*(x@z);width=np.linalg.norm(x)
        if width<1e-5:raise ValueError('The two source courses must be separated in their reference plane')
        x/=width;return np.column_stack([x,np.cross(z,x),z]),width
    before,old_width=frame(a,source_axis);after,new_width=frame(b,target_axis);width_scale=new_width/old_width
    if not .5<=width_scale<=2.5 or not .5<=axial_scale<=2.5:raise ValueError('Explicit vascular carriage scales outside [.5,2.5]')
    R=after@before.T;A=after@np.diag([width_scale,width_scale,axial_scale])@before.T
    old,new=a.mean(0),b.mean(0);q=new+(v-old)@A.T;physical=np.asarray(physical_ids,int);owned=np.zeros(len(v),bool)
    if physical.shape!=(len(v),):raise ValueError('Original physical identities required')
    for section in sections:
        s=np.asarray(section,int)
        if owned[s].any():raise ValueError('Original section ownership overlaps')
        centre=v[s].mean(0);point=new+A@(centre-old)
        for raw in s:
            aliases=physical==physical[raw];q[aliases]=(v[aliases]-centre)@R.T+point;owned[aliases]=True
    return q,np.repeat(R[None],len(v),axis=0),dict(method='Paired palmar-vessel reference frame with positive longitudinal and transverse carriage',
        sourcePair_m=a.tolist(),targetPair_m=b.tolist(),sourceAxis=np.asarray(source_axis).tolist(),targetAxis=np.asarray(target_axis).tolist(),
        transverseCarriageScale=float(width_scale),axialCarriageScale=float(axial_scale),sourceSectionCalibreUnscaled=True,
        carrierDeterminant=float(np.linalg.det(A)),sourceSkeletonChanged=False,targetSkeletonChanged=False,anatomicalAcceptance=False)


def surface_frame_carriage(vertices,sections,carriers,physical_ids,*,maximum_carrier_scale=2.5):
    """Use full named-bone surface frames, preserving each source section.

    Affine width/length ratios position centres only. Proper blended rotations
    carry their unchanged source offsets. Source bones stay read-only guides.
    """
    from scipy.spatial.transform import Rotation
    v=np.asarray(vertices,float);physical=np.asarray(physical_ids,int)
    if not carriers or physical.shape!=(len(v),):raise ValueError('Explicit named-bone carriers and physical identities required')
    if not 1<=maximum_carrier_scale<=4:raise ValueError('Explicit finite carrier scale bound within [1,4] required')
    transforms=[];quaternions=[]
    for c in carriers:
        a,b=np.asarray(c['sourceFrame']),np.asarray(c['targetFrame']);scale=np.asarray(c['scales'])
        if scale.shape!=(3,)or not np.isfinite(scale).all()or np.any(scale<.5)or np.any(scale>maximum_carrier_scale):
            raise ValueError('Measured finite bone carrier dimensions outside the declared bound')
        for f in [a,b]:
            if not np.allclose(f.T@f,np.eye(3),atol=1e-11,rtol=0)or np.linalg.det(f)<1-1e-11:
                raise ValueError('Proper source and recipient bone surface frames required')
        transforms.append(b@np.diag(scale)@a.T);quaternions.append(Rotation.from_matrix(b@a.T).as_quat())
    quaternions=np.array(quaternions)
    quaternions[np.sum(quaternions*quaternions[0],axis=1)<0]*=-1
    def mapped(points):
        distance=np.column_stack([c['index'].closest_points(points)[0]['distance_m']for c in carriers])
        weights=1/(distance+.002)**4;weights/=weights.sum(1)[:,None]
        placed=np.stack([np.asarray(c['targetCentre'])+(points-c['sourceCentre'])@M.T for c,M in zip(carriers,transforms)],axis=1)
        q=np.einsum('nk,nki->ni',weights,placed)
        quaternion=weights@quaternions;quaternion/=np.linalg.norm(quaternion,axis=1)[:,None]
        return q,Rotation.from_quat(quaternion).as_matrix()
    q,R=mapped(v);owned=np.zeros(len(v),bool)
    # Query corresponding section centres as one batch rather than making
    # one native spatial query per vessel ring and named bone.
    centres=np.array([v[np.asarray(s,int)].mean(0)for s in sections]).reshape(-1,3)
    placed,rotations=mapped(centres)if len(centres)else(np.empty((0,3)),np.empty((0,3,3)))
    for section,centre,point,rotation in zip(sections,centres,placed,rotations):
        s=np.asarray(section,int)
        if owned[s].any():raise ValueError('Source section ownership overlaps')
        for raw in s:
            aliases=physical==physical[raw];q[aliases]=(v[aliases]-centre)@rotation.T+point;R[aliases]=rotation;owned[aliases]=True
    return q,R,dict(method='Named-bone surface-frame partition of unity; proper source section poses',
                    namedBoneCarriers=len(carriers),distanceWeightRegularizer_m=.002,distanceWeightExponent=4,
                    measuredScales=[np.asarray(c['scales']).tolist()for c in carriers],sourceSectionCalibreUnscaled=True,
                    maximumCarrierScale=maximum_carrier_scale,
                    sourceSkeletonChanged=False,targetSkeletonChanged=False,anatomicalAcceptance=False)


def paired_bone_transport(vertices,sections,source_bones,target_bones,source_palm_normal,target_palm_normal,*,physical_ids=None,transverse_station_scales=None):
    v=np.asarray(vertices,float);a=np.asarray(source_bones,float);b=np.asarray(target_bones,float)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or len(a)<2:
        raise ValueError('Explicit matching named bone stations required')
    if not all(np.isfinite(x).all()for x in [v,a,b]):raise ValueError('Finite source coordinates required')
    if np.any(np.linalg.norm(np.diff(a,axis=0),axis=1)<1e-8)or np.any(np.linalg.norm(np.diff(b,axis=0),axis=1)<1e-8):
        raise ValueError('Distinct corresponding bone stations required')
    axes=[];frames=[]
    t=np.linspace(-.75,len(a)-.25,257)
    for bones,normal in [(a,source_palm_normal),(b,target_palm_normal)]:
        controls=np.r_[(2*bones[0]-bones[1])[None],bones,(2*bones[-1]-bones[-2])[None]]
        curve=PchipInterpolator(np.arange(-1,len(bones)+1),controls)
        axis=curve(t);tangent=unit(curve.derivative()(t))
        axes.append(axis);frames.append(rotation_minimizing_frames(axis,tangent,initial_reference=normal))
    rotations=frames[1]@frames[0].transpose(0,2,1)
    scales=None
    if transverse_station_scales is not None:
        scales=np.asarray(transverse_station_scales,float)
        if scales.shape!=(len(a),2)or not np.isfinite(scales).all()or np.any(scales<.5)or np.any(scales>2.5):
            raise ValueError('Explicit paired surface-width ratios within [.5,2.5] required')

    def map_points(points):
        q,R=axis_transport(points,*axes,rotations)
        if scales is not None:
            from scipy.spatial.transform import Rotation,Slerp
            axis=axes[0];segments=np.diff(axis,axis=0)
            fraction=np.clip(np.einsum('pki,ki->pk',points[:,None]-axis[:-1],segments)/np.sum(segments*segments,axis=1),0,1)
            closest=axis[:-1][None]+fraction[:,:,None]*segments[None]
            k=np.argmin(np.sum((points[:,None]-closest)**2,axis=2),axis=1);u=fraction[np.arange(len(points)),k]
            station=np.interp(k+u,np.arange(len(t)),t)
            ratio=np.column_stack([np.interp(station,np.arange(len(a)),scales[:,j])for j in range(2)])
            before=Slerp(np.arange(len(t)),Rotation.from_matrix(frames[0]))(k+u).as_matrix()
            offset=points-closest[np.arange(len(points)),k]
            local=np.einsum('nji,nj->ni',before,offset)
            local[:,:2]*=ratio-1.;local[:,2]=0.
            q+=np.einsum('nij,nj->ni',R,np.einsum('nij,nj->ni',before,local))
        return q,R

    q,R=map_points(v)
    physical=np.asarray(physical_ids,int)if physical_ids is not None else np.arange(len(v))
    if physical.shape!=(len(v),):raise ValueError('One source physical identity per raw vertex required')
    owned=np.zeros(len(v),bool)
    for section in sections:
        s=np.asarray(section,int)
        if owned[s].any():raise ValueError('Rigid original section ownership overlaps')
        centre=v[s].mean(0);point,rotation=map_points(centre[None])
        for raw in s:
            aliases=physical==physical[raw]
            q[aliases]=(v[aliases]-centre)@rotation[0].T+point[0];R[aliases]=rotation[0];owned[aliases]=True
    if np.array_equal(a,b)and np.array_equal(source_palm_normal,target_palm_normal)and(scales is None or np.all(scales==1.)):
        q=v.copy();R=np.broadcast_to(np.eye(3),(len(v),3,3)).copy()
    return q,R,dict(namedSourceStations=a.tolist(),namedTargetStations=b.tolist(),
                    rigidSections=len(sections),fixedSkeleton=True,anatomicalAcceptance=False,
                    measuredTransverseStationScales=None if scales is None else scales.tolist(),
                    widthScalingActsOnSectionCentresOnly=True,
                    mapping='Corresponding bone station parameters and proper transported section frames')
