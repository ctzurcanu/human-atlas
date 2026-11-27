"""Recover occupied and fluid envelopes of a finite, annular-ended tube mesh.

No anatomical dimensions are invented. Native side faces and annular caps
remain literal. Only temporary fluid-mouth disks are added. The same disks
belong to both envelope caps and cancel when their difference is formed.
"""
from dataclasses import dataclass
import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


class PipeEnvelopeError(ValueError):
    def __init__(self,message,diagnostics=None):
        super().__init__(message);self.diagnostics=diagnostics


@dataclass
class PipeEnvelopes:
    positions:np.ndarray
    original_faces:np.ndarray
    outer_faces:np.ndarray
    fluid_faces:np.ndarray
    cap_face_ids:np.ndarray
    inner_side_face_ids:np.ndarray
    outer_side_face_ids:np.ndarray
    fluid_disk_faces:np.ndarray
    report:dict


def face_adjacency(faces):
    f=np.asarray(faces,int);edges=np.vstack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);owners=np.tile(np.arange(len(f)),3)
    _,inverse,count=np.unique(np.sort(edges,axis=1),axis=0,return_inverse=True,return_counts=True)
    order=np.argsort(inverse,kind='stable');starts=np.r_[0,np.cumsum(count)[:-1]];paired=starts[count==2]
    return np.column_stack([owners[order[paired]],owners[order[paired+1]]]),count


def face_components(faces,adjacency=None):
    f=np.asarray(faces,int);pairs=face_adjacency(f)[0]if adjacency is None else adjacency
    graph=coo_matrix((np.ones(2*len(pairs)),(np.r_[pairs[:,0],pairs[:,1]],np.r_[pairs[:,1],pairs[:,0]])),shape=(len(f),len(f))).tocsr()
    return connected_components(graph,directed=False)


def boundary_loops(faces):
    f=np.asarray(faces,int);e=np.vstack([f[:,[0,1]],f[:,[1,2]],f[:,[2,0]]]);_,inv,count=np.unique(np.sort(e,axis=1),axis=0,return_inverse=True,return_counts=True);b=e[count[inv]==1]
    successor=dict(map(tuple,b.tolist()))
    if len(successor)!=len(b)or set(successor)!=set(successor.values()):raise PipeEnvelopeError('Boundary is not a disjoint set of directed cycles')
    unseen=set(successor);loops=[]
    while unseen:
        loop=[min(unseen)]
        while successor[loop[-1]]!=loop[0]:
            loop.append(successor[loop[-1]])
            if len(loop)>len(b):raise PipeEnvelopeError('Boundary cycle did not close')
        unseen.difference_update(loop);loops.append(np.array(loop,int))
    return loops


def projected_area(points):
    p=np.asarray(points,float);return float(np.linalg.norm(np.cross(p-p.mean(0),np.roll(p,-1,axis=0)-p.mean(0)).sum(0))/2)


def signed_volume(v,f):
    p=np.asarray(v,float)[np.asarray(f,int)];p=p-p.mean(axis=(0,1));return float(np.einsum('ti,ti->t',p[:,0],np.cross(p[:,1],p[:,2])).sum()/6)


def close_with_shared_disks(v,f,expected_mouths=2):
    """Close a declared number of mouths with computational fluid faces.

    Ordinary pipe callers retain the strict two-mouth default. A branched
    reservoir must explicitly declare its different source opening count.
    """
    f=np.asarray(f)
    if not np.issubdtype(f.dtype,np.integer)or(f.size and (f.min()<0 or f.max()>np.iinfo(np.int64).max)):raise PipeEnvelopeError('Finite integer face indices required')
    f=f.astype(np.int64,copy=False)
    loops=boundary_loops(f)
    if len(loops)!=expected_mouths:raise PipeEnvelopeError(f'Side surface must have exactly {expected_mouths} declared mouths')
    points=list(v);disks=[]
    for loop in loops:
        center=len(points);points.append(v[loop].mean(0))
        disks.extend([[int(b),int(a),center]for a,b in zip(loop,np.roll(loop,-1))])
    disks=np.asarray(disks,int);return np.asarray(points),np.vstack([f,disks]),disks


def recover_paired_sheets(v,f,sheet_labels,*,expected_mouths=2):
    """Use explicit source shell lineage, including non-planar annular ends.

    Labels must come from source topology, never guessed spatial proximity.
    Both side sheets and both annular patches are independently validated.
    The caller must also verify fluid containment and wall reconstruction.
    """
    if not isinstance(expected_mouths,int)or not 2<=expected_mouths<=512:
        raise PipeEnvelopeError('An explicit bounded source mouth count >= 2 is required')
    v=np.asarray(v,float);f=np.asarray(f,int);labels=np.asarray(sheet_labels,int)
    if labels.shape!=(len(v),)or set(labels)!={0,1}:raise PipeEnvelopeError('Two explicit source sheet labels are required')
    v,first,inv=np.unique(v,axis=0,return_index=True,return_inverse=True)
    if not np.array_equal(labels[first][inv],labels):raise PipeEnvelopeError('Opposite sheet labels share a collapsed coordinate')
    labels=labels[first];f=inv[f];pairs,count=face_adjacency(f)
    if np.any(count!=2)or face_components(f,pairs)[0]!=1:raise PipeEnvelopeError('Source is not one closed manifold wall')
    mixed=np.any(labels[f]!=labels[f][:,:1],axis=1);capids=np.flatnonzero(mixed);cap=f[capids]
    capcount,caplabels=face_components(cap)
    if capcount!=expected_mouths:raise PipeEnvelopeError('Source lineage must expose the declared number of separate annular ends')
    caploops=[]
    for k in range(capcount):
        ff=cap[caplabels==k];loops=boundary_loops(ff)
        if len(loops)!=2:raise PipeEnvelopeError('Each source end must have two annular boundary cycles')
        if not all(len(set(labels[loop]))==1 for loop in loops)or labels[loops[0][0]]==labels[loops[1][0]]:raise PipeEnvelopeError('Annular boundaries must belong to opposite source sheets')
        caploops.extend(loops)
    sides=[]
    for label in [0,1]:
        ids=np.flatnonzero((~mixed)&(labels[f[:,0]]==label));side=f[ids]
        if face_components(side)[0]!=1:raise PipeEnvelopeError('Each source side must be connected')
        loops=boundary_loops(side)
        if len(loops)!=expected_mouths or {frozenset(loop)for loop in loops}!={frozenset(loop)for loop in caploops if labels[loop[0]]==label}:raise PipeEnvelopeError('Side mouths and annular ends do not match exactly')
        xyz,closed,_=close_with_shared_disks(v,side,expected_mouths=expected_mouths);sides.append((ids,signed_volume(xyz,closed)))
    inner=[s for s in sides if s[1]<0];outer=[s for s in sides if s[1]>0]
    if len(inner)!=1 or len(outer)!=1:raise PipeEnvelopeError('Source sides do not have opposite inner/outer orientation')
    innerids,outerids=inner[0][0],outer[0][0]
    positions,fluid,disks=close_with_shared_disks(v,f[innerids][:,[0,2,1]],expected_mouths=expected_mouths)
    exterior=np.vstack([f[np.setdiff1d(np.arange(len(f)),innerids)],disks])
    for name,faces in [('fluid',fluid),('outer',exterior)]:
        if np.any(face_adjacency(faces)[1]!=2)or signed_volume(positions,faces)<=0:raise PipeEnvelopeError(name+' envelope does not close consistently')
    report=dict(method='explicit paired source-sheet topology',sourceVertices=len(v),sourceFaces=len(f),annularCapFaces=len(capids),mouths=expected_mouths,innerSideFaces=len(innerids),outerSideFaces=len(outerids),temporaryFluidDiskTriangles=len(disks),
        nativeCoordinatesAndSideFacesLiteral=True,annularTissueCapsLiteral=True,temporaryFluidDisksSharedExactly=True,nonPlanarNativeAnnularEndsRetained=True,anatomicallyAccepted=False)
    return PipeEnvelopes(positions,f,exterior,fluid,capids,innerids,outerids,disks,report)


def fluid_through_mouths(envelopes):
    """Extend only computational fluid mouths through their native end rims.

    The cutter ends beyond the entire annular rim, avoiding coincident end
    disks and retaining a through opening after outer-minus-fluid booleans.
    No native outer-envelope point or inner-side face moves.
    """
    r=envelopes;v=r.positions;inner=r.original_faces[r.inner_side_face_ids][:,[0,2,1]];loops=boundary_loops(inner)
    caps=r.original_faces[r.cap_face_ids];nc,labels=face_components(caps);assert nc==2
    points=list(v);faces=list(inner);rows=[]
    for loop in loops:
        cap=next(caps[labels==i]for i in range(nc)if set(loop)<=set(caps[labels==i].ravel()))
        caploops=boundary_loops(cap);best=min(caploops,key=lambda p:np.linalg.svd(v[p]-v[p].mean(0),compute_uv=False)[-1])
        p=v[best];_,_,vectors=np.linalg.svd(p-p.mean(0),full_matrices=False);direction=vectors[-1]
        touching=inner[np.any(np.isin(inner,loop),axis=1)];inside=np.setdiff1d(np.unique(touching),loop)
        if not len(inside):raise PipeEnvelopeError('Cannot establish the inward side of a mouth')
        direction*=(-1 if np.dot(direction,v[inside].mean(0)-v[loop].mean(0))>0 else 1)
        center=v[loop].mean(0);cap_points=v[np.unique(cap)]
        guard=2*float(np.linalg.norm(np.spacing(cap_points.astype('<f4')).astype(float),axis=1).max())
        distance=max(0.,float(((cap_points-center)@direction).max()-((v[loop]-center)@direction).min()))+guard
        new=np.arange(len(points),len(points)+len(loop));points.extend(v[loop]+distance*direction);tip=len(points);points.append(v[loop].mean(0)+distance*direction)
        for i,j in zip(range(len(loop)),np.roll(np.arange(len(loop)),-1)):
            a,b,aa,bb=int(loop[i]),int(loop[j]),int(new[i]),int(new[j]);faces.extend([[b,a,aa],[b,aa,bb],[bb,aa,tip]])
        rows.append(dict(extensionMeters=distance,numericalGuardMeters=guard,direction=direction.tolist(),entireEndRimBehindCutterCap=True))
    xyz,ff=np.asarray(points),np.asarray(faces,int)
    if np.any(face_adjacency(ff)[1]!=2)or signed_volume(xyz,ff)<=0:raise PipeEnvelopeError('Extended fluid cutter does not close consistently')
    return xyz,ff,rows


def recover(v,f,*,angular_tolerance_rad=.02,planarity_tolerance=None,expected_mouths=2):
    """One closed hollow component with explicitly declared annular mouths.

    angular_tolerance_rad groups candidate patches only; it cannot waive the
    separate plane residual, annulus topology or envelope-difference checks.
    Default planarity tolerance is derived from two stored float32 spacings.
    """
    if not isinstance(expected_mouths,int)or not 2<=expected_mouths<=512:
        raise PipeEnvelopeError('An explicit bounded source mouth count >= 2 is required')
    v=np.asarray(v,float);f=np.asarray(f,int);v,inv=np.unique(v,axis=0,return_inverse=True);f=inv[f]
    pairs,count=face_adjacency(f)
    if np.any(count!=2):raise PipeEnvelopeError('Source must be a closed two-manifold')
    ncomp,_=face_components(f,pairs)
    if ncomp!=1:raise PipeEnvelopeError('Recover one connected source component at a time')
    tri=v[f];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);length=np.linalg.norm(normal,axis=1)
    if np.any(length==0):raise PipeEnvelopeError('Degenerate source face')
    normal/=length[:,None]
    tol=float(2*np.linalg.norm(np.spacing(v.astype('<f4')).astype(float),axis=1).max())if planarity_tolerance is None else float(planarity_tolerance)
    if not np.isfinite(tol)or tol<=0:raise PipeEnvelopeError('Positive finite plane tolerance required')
    # Thin annular facets amplify position rounding into angular uncertainty.
    # Bound that amplification from each stored corner's float32 spacing;
    # the independent whole-patch plane residual remains mandatory.
    point_error=np.linalg.norm(np.spacing(tri.astype('<f4')).astype(float),axis=2)/2
    edge1,edge2=tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]
    error1,error2=point_error[:,0]+point_error[:,1],point_error[:,0]+point_error[:,2]
    cross_error=np.linalg.norm(edge1,axis=1)*error2+np.linalg.norm(edge2,axis=1)*error1+error1*error2
    angular_error=np.arcsin(np.minimum(1.,cross_error/length))
    angle_limit=np.minimum(np.pi,angular_tolerance_rad+angular_error[pairs[:,0]]+angular_error[pairs[:,1]])
    same=np.einsum('ij,ij->i',normal[pairs[:,0]],normal[pairs[:,1]])>=np.cos(angle_limit)
    npatch,labels=face_components(f,pairs[same]);caps=[];candidates=[]
    for patch in range(npatch):
        ids=np.flatnonzero(labels==patch)
        if len(ids)<6:continue
        ff=f[ids];loops=boundary_loops(ff)
        if len(loops)!=2:continue
        used=np.unique(ff);edgecount=len(np.unique(np.sort(np.vstack([ff[:,[0,1]],ff[:,[1,2]],ff[:,[2,0]]]),axis=1),axis=0))
        if len(used)-edgecount+len(ff)!=0:continue
        points=v[used];_,_,vectors=np.linalg.svd(points-points.mean(0),full_matrices=False);residual=float(abs((points-points.mean(0))@vectors[-1]).max())
        candidates.append(dict(patch=int(patch),faces=len(ids),loops=[len(loop)for loop in loops],planeResidual=residual,planeTolerance=tol,center=points.mean(0).tolist(),normal=vectors[-1].tolist()))
        if residual>tol:continue
        loops.sort(key=lambda loop:projected_area(v[loop]))
        caps.append(dict(ids=ids,inner=loops[0],outer=loops[1],planeResidual=residual))
    if len(caps)!=expected_mouths:raise PipeEnvelopeError(f'Expected {expected_mouths} planar annular end patches; found {len(caps)}',dict(annularCandidates=candidates,patches=int(npatch)))
    capids=np.concatenate([c['ids']for c in caps]);sideids=np.setdiff1d(np.arange(len(f)),capids);side=f[sideids]
    ns,sl=face_components(side)
    if ns!=2:raise PipeEnvelopeError(f'Expected two separate side surfaces; found {ns}')
    innerids=[];outerids=[]
    for label in range(ns):
        ids=sideids[sl==label];used=set(f[ids].ravel());is_inner=all(set(c['inner'])<=used for c in caps);is_outer=all(set(c['outer'])<=used for c in caps)
        if is_inner==is_outer:raise PipeEnvelopeError('Side surfaces do not match both annular end boundaries')
        (innerids if is_inner else outerids).extend(ids)
    innerids=np.asarray(innerids,int);outerids=np.asarray(outerids,int)
    inner=f[innerids][:,[0,2,1]];loops=boundary_loops(inner)
    if len(loops)!=expected_mouths:raise PipeEnvelopeError('Inner surface must have the declared number of mouths')
    points=list(v);disks=[]
    for loop in loops:
        center=len(points);points.append(v[loop].mean(0))
        disks.extend([[int(b),int(a),center]for a,b in zip(loop,np.roll(loop,-1))])
    positions=np.asarray(points);disks=np.asarray(disks,int);fluid=np.vstack([inner,disks]);outer=np.vstack([f[np.setdiff1d(np.arange(len(f)),innerids)],disks])
    for name,faces in [('fluid',fluid),('outer',outer)]:
        if np.any(face_adjacency(faces)[1]!=2):raise PipeEnvelopeError(name+' envelope is not closed')
        if signed_volume(positions,faces)<=0:raise PipeEnvelopeError(name+' envelope orientation is invalid')
    report=dict(sourceVertices=len(v),sourceFaces=len(f),annularCapFaces=len(capids),mouths=expected_mouths,innerSideFaces=len(innerids),outerSideFaces=len(outerids),temporaryFluidDiskTriangles=len(disks),
        maximumCapPlaneResidual=float(max(c['planeResidual']for c in caps)),planeTolerance=tol,angularPatchGroupingToleranceRadians=float(angular_tolerance_rad),maximumFloat32NormalAngleUncertaintyRadians=float(angular_error.max()),
        nativeCoordinatesAndSideFacesLiteral=True,annularTissueCapsLiteral=True,temporaryFluidDisksSharedExactly=True,anatomicallyAccepted=False)
    return PipeEnvelopes(positions,f,outer,fluid,capids,innerids,outerids,disks,report)
