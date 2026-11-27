"""Cancel certified internal opposed facets in an otherwise regular open tube.

Actual source openings stay open. Existing-ring disks close only a separate
computational mesh for winding queries; none enters the returned face bank.
Queries must lie outside every ring's enclosing ball, so alternative fans
inside those same convex hulls cannot change the winding classification.
"""
from collections import Counter
import time
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.contact import pair_contacts
from geometry_tools.spatial import BVH
from .tube_continuity import boundary_loops
from .convex_surface_channel import earcut


def _boundary_chain(faces):
    counts=Counter()
    for tri in faces:
        for a,b in zip(tri,np.roll(tri,-1)):
            a,b=int(a),int(b);counts[min(a,b),max(a,b)]+=1 if a<b else -1
    return {k:v for k,v in counts.items() if v}


def cancel_open_envelope_internal_facets(mesh,*,wall_seconds=6.):
    import igl
    start=time.monotonic()
    def deadline():
        if time.monotonic()-start>wall_seconds:
            raise TimeoutError('Open-source internal-facet proof deadline')
    physical,alias=np.unique(mesh.vertices,axis=0,return_inverse=True)
    f=alias[mesh.faces]
    _,group,count=np.unique(np.sort(f,axis=1),axis=0,return_inverse=True,return_counts=True)
    pairs=[]
    for k in np.flatnonzero(count==2):
        ids=np.flatnonzero(group==k);a,b=ids
        if any(np.array_equal(f[a],np.roll(f[b][::-1],j)) for j in range(3)):
            pairs.append(ids.tolist())
    if not pairs:raise ValueError('No exactly opposed paired source facets')
    removed=np.asarray(pairs,int).ravel();keep=np.setdiff1d(np.arange(len(f)),removed)
    cf=f[keep]
    candidate=Mesh(mesh.vertices,mesh.faces[keep],mesh.frame,mesh.attributes,mesh.source)
    d=diagnostics(candidate,include_pairs=True);deadline()
    if d['nonmanifold_edges'] or d['nonmanifold_vertices'] or d['self_intersection_pairs'] or not d['orientable']:
        raise ValueError('Retained open envelope must be manifold, orientable and self-disjoint')
    if _boundary_chain(f)!=_boundary_chain(cf):
        raise ValueError('Source oriented boundary chain changed')
    loops=boundary_loops(cf)
    if not loops:raise ValueError('Use the closed-envelope proof for a closed source')
    caps=[];balls=[]
    for loop in loops:
        deadline();ring=np.asarray(loop,int)[::-1];p=physical[ring];center=p.mean(0)
        _,_,basis=np.linalg.svd(p-center,full_matrices=False);xy=(p-center)@basis[:2].T
        coordinates={int(i):q for i,q in zip(ring,xy)}
        fan=np.asarray(earcut(ring,coordinates),int)
        if fan.shape!=(len(ring)-2,3):
            raise ValueError('Computational boundary triangulation lost a source control')
        tri=physical[fan]
        if np.any(np.linalg.norm(np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]),axis=1)<=1e-20):
            raise ValueError('Degenerate computational boundary fan')
        caps.append(fan);balls.append((center,float(np.linalg.norm(p-center,axis=1).max())))
    cap_faces=np.concatenate(caps);closed_faces=np.r_[cf,cap_faces]
    closed=Mesh(physical,closed_faces,mesh.frame)
    dc=diagnostics(closed,include_pairs=True);deadline()
    if dc['nonmanifold_edges'] or dc['nonmanifold_vertices'] or dc['self_intersection_pairs'] or not dc['orientable'] or boundary_loops(closed_faces):
        raise ValueError('Computational closure must be closed, manifold and self-disjoint')
    tris=mesh.triangles[np.asarray(pairs)[:,0]]
    normals=np.cross(tris[:,1]-tris[:,0],tris[:,2]-tris[:,0]);length=np.linalg.norm(normals,axis=1)
    if np.any(length<=1e-20):raise ValueError('Degenerate opposed source facet')
    normals/=length[:,None];centers=tris.mean(1)
    distance=BVH(candidate).closest_points(centers)[0]['distance_m']
    if np.any(distance<1e-7):raise ValueError('Removed facet touches the retained source exterior')
    epsilon=np.minimum(distance*.1,1e-5)
    points=np.concatenate([centers,centers+normals*epsilon[:,None],centers-normals*epsilon[:,None]])
    margins=[float((np.linalg.norm(points-center,axis=1)-radius).min()) for center,radius in balls]
    if min(margins)<=1e-7:
        raise ValueError('Winding probes overlap a source-opening hull; closure-independent interior is unproven')
    origin=physical.mean(0)
    w=igl.winding_number(physical-origin,closed_faces,points-origin)
    if np.max(np.abs(np.abs(w)-1))>1e-9:
        raise ValueError('Opposed source facets are not inside one retained vascular envelope')
    # Add the same computational fans to the old oriented chain. Exact source
    # facet cancellation must preserve winding on both sides of each facet.
    old_w=igl.winding_number(physical-origin,np.r_[f,cap_faces],points[len(centers):]-origin)
    if np.max(np.abs(w[len(centers):]-old_w))>1e-9:
        raise ValueError('Source winding changed')
    deleted=Mesh(mesh.vertices,mesh.faces[removed],mesh.frame)
    _,contact=pair_contacts(deleted,candidate,0.,max_tests=250000);deadline()
    if not contact['complete_within_pair_budget'] or contact['crossing_face_pairs'] or contact['coplanar_overlap_face_pairs']:
        raise ValueError('Removed facets intersect the retained source exterior')
    _,cap_contact=pair_contacts(deleted,Mesh(physical,cap_faces,mesh.frame),0.,max_tests=250000);deadline()
    if not cap_contact['complete_within_pair_budget'] or cap_contact['crossing_face_pairs'] or cap_contact['coplanar_overlap_face_pairs']:
        raise ValueError('Removed facets intersect a computational mouth closure')
    return keep,dict(complete=True,removedSourceFaceIds=removed.tolist(),retainedSourceFaceIds=keep.tolist(),
        opposedSourcePairs=pairs,sourceBoundaryChainLiteral=True,sourceCoordinatesAndFieldsLiteral=True,
        retainedEnvelopeManifold=True,retainedEnvelopeSelfPairs=0,sourceOpenings=len(loops),
        computationalCaps=len(cap_faces),computationalCapFaces=cap_faces.tolist(),computationCapsEmitted=0,
        computationalDiskMethod='Existing source-boundary ear clipping; no inserted controls',
        minimumRemovedCentroidDistanceToExterior_m=float(distance.min()),
        minimumProbeClearanceFromAllMouthHullBalls_m=min(margins),
        maximumInteriorWindingError=float(np.abs(np.abs(w)-1).max()),
        maximumSourceWindingDifference=float(np.abs(w[len(centers):]-old_w).max()),
        allOtherTriangleContactsCanOnlyBeRemoved=True,finiteWallOrClinicalPatencyEstablished=False,
        anatomicalAcceptance=False)
