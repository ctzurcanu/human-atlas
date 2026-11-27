"""One local arteriovenous carrier around explicitly corresponding digit bones.

The neighbouring source digits bound its support. Source sections keep their
shape under proper rotations; the recipient skeleton is an immutable input.
This is a placement proposal, not tissue-clearance or anatomical acceptance.
"""
import numpy as np
from scipy.interpolate import PchipInterpolator
from scipy.spatial.transform import Rotation,Slerp
from .tube_sections import rotation_minimizing_frames,unit
from .tube_axis_map import axis_transport


def polyline_distance(points,curve):
    p=np.asarray(points,float);c=np.asarray(curve,float);d=np.diff(c,axis=0)
    length=np.sum(d*d,axis=1)
    if np.any(length<1e-16):raise ValueError('Distinct bone stations required')
    u=np.clip(np.einsum('nki,ki->nk',p[:,None]-c[:-1],d)/length,0,1)
    nearest=c[:-1][None]+u[:,:,None]*d[None]
    return np.linalg.norm(p[:,None]-nearest,axis=2).min(1)


def local_bone_cohort(vertices,sections,source_bones,target_bones,
                      source_palm_normal,target_palm_normal,neighbour_bones,*,
                      palm_blend=(-.4,.6),separation_blend_m=.008,radial_scale=None,
                      baseline_vertices=None,neighbour_partition=True,pose_support=None):
    v=np.asarray(vertices,float);a=np.asarray(source_bones,float);b=np.asarray(target_bones,float)
    if a.shape!=b.shape or a.ndim!=2 or a.shape[1]!=3 or len(a)<3:
        raise ValueError('Three or more explicit corresponding bone stations required')
    if not all(np.isfinite(x).all()for x in [v,a,b])or not neighbour_bones:
        raise ValueError('Finite source geometry and explicit neighbouring digits required')
    baseline=v if baseline_vertices is None else np.asarray(baseline_vertices,float)
    if baseline.shape!=v.shape or not np.isfinite(baseline).all():
        raise ValueError('Explicit corresponding installed controls required')
    lo,hi=map(float,palm_blend)
    if not -1<=lo<hi<=1 or not .001<=separation_blend_m<=.02:
        raise ValueError('Bounded explicit support transition required')
    t=np.linspace(-.75,len(a)-.25,257);axes=[];frames=[]
    for bones,normal in [(a,source_palm_normal),(b,target_palm_normal)]:
        controls=np.r_[(2*bones[0]-bones[1])[None],bones,(2*bones[-1]-bones[-2])[None]]
        curve=PchipInterpolator(np.arange(-1,len(bones)+1),controls)
        axis=curve(t);axes.append(axis)
        frames.append(rotation_minimizing_frames(axis,unit(curve.derivative()(t)),initial_reference=normal))
    rotations=frames[1]@frames[0].transpose(0,2,1)
    def map_points(p,base):
        q,R=axis_transport(p,axes[0],axes[1],rotations)
        d=np.diff(axes[0],axis=0);u=np.clip(np.einsum('nki,ki->nk',p[:,None]-axes[0][:-1],d)/np.sum(d*d,axis=1),0,1)
        near=axes[0][:-1][None]+u[:,:,None]*d[None]
        k=np.argmin(np.sum((p[:,None]-near)**2,axis=2),axis=1)
        station=np.interp(k+u[np.arange(len(p)),k],np.arange(len(t)),t)
        radial=np.linalg.norm(p-near[np.arange(len(p)),k],axis=1)
        other=np.column_stack([polyline_distance(p,c)for c in neighbour_bones]).min(1)
        w=np.clip((other-radial)/separation_blend_m,0,1)if neighbour_partition else np.ones(len(p))
        longitudinal=np.clip((station-lo)/(hi-lo),0,1)
        smooth=lambda x:x*x*(3-2*x)
        w=smooth(w)*smooth(longitudinal)
        if pose_support is not None:
            support=np.asarray(pose_support(p),bool)
            if support.shape!=(len(p),):raise ValueError('One explicit pose-domain flag per query required')
            w[~support]=0.
        if radial_scale is not None:
            live=w>0;delta=p[live]-near[np.arange(len(p)),k][live]
            scale=np.asarray(radial_scale(station[live],delta),float)
            if scale.shape!=(int(live.sum()),)or not np.isfinite(scale).all()or np.any(scale<=0):
                raise ValueError('Finite positive radial placement factors required')
            q[live]+=np.einsum('nij,nj->ni',R[live],delta)*(scale-1)[:,None]
        placed=base+w[:,None]*(q-base)
        rotation=Rotation.from_rotvec(Rotation.from_matrix(R).as_rotvec()*w[:,None]).as_matrix()
        placed[w==0]=base[w==0];rotation[w==0]=np.eye(3)
        return placed,rotation,w
    q,R,w=map_points(v,baseline);owned=np.zeros(len(v),bool);centres=np.array([v[s].mean(0)for s in sections])
    if len(sections):
        base_centres=np.array([baseline[s].mean(0)for s in sections])
        target,rotations,weights=map_points(centres,base_centres)
        for s,c,p,rot,weight in zip(sections,centres,target,rotations,weights):
            if owned[s].any():raise ValueError('Explicit disjoint physical source sections required')
            if baseline_vertices is not None:
                u,_,vh=np.linalg.svd((v[s]-c).T@(baseline[s]-baseline[s].mean(0)))
                sign=np.eye(3);sign[-1,-1]=np.linalg.det(vh.T@u.T)
                base_rot=vh.T@sign@u.T
                if weight==0:
                    q[s]=baseline[s];R[s]=base_rot;owned[s]=True;w[s]=0.;continue
                full_rot=Rotation.from_rotvec(Rotation.from_matrix(rot).as_rotvec()/weight).as_matrix()
                rot=Slerp([0.,1.],Rotation.from_matrix(np.stack([base_rot,full_rot])))([weight]).as_matrix()[0]
            q[s]=(v[s]-c)@rot.T+p;R[s]=rot;w[s]=weight;owned[s]=True
    return q,R,dict(method='Corresponding digit bone station frames with neighbouring-digit support and proper section rotations',
                    sourceSectionCalibreUnchanged=baseline_vertices is None,sourceBranchesRemoved=0,recipientSkeletonChanged=False,
                    palmBlendStations=[lo,hi],separationBlend_m=separation_blend_m,
                    nonzeroControls=int((w>0).sum()),fullWeightControls=int((w>=1-1e-12).sum()),
                    maximumControlMotion_m=float(np.linalg.norm(q-v,axis=1).max(initial=0.)),
                    measuredBoneRadialPlacement=radial_scale is not None,
                    installedBaselineOutsideSupportPreserved=baseline_vertices is not None,
                    independentNeighbourPartition=neighbour_partition,
                    anatomicalAcceptance=False,fullTriangleAndTissueAuditRequired=True)
