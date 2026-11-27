"""Remove explicitly classified interior source partitions, preserving exteriors.

This cannot classify a valve or other real tissue as an artifact. Callers must
declare the source facets, establish their ownership and review actual cuts.
Computational mouth closures remain separate and are never emitted.
"""
import time
import numpy as np
from geometry_tools.model import Mesh
from geometry_tools.open3d_backend import diagnostics
from geometry_tools.spatial import BVH
from geometry_tools.contact import pair_contacts
from .tube_continuity import boundary_loops
from .convex_surface_channel import earcut
from .shared_boundary_arc import topology

def edges(f):
    return np.unique(np.sort(f[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1),axis=0,return_counts=True)

def remove_internal_source_sheet(mesh,face_ids,*,wall_seconds=6.):
    import igl
    tick=time.monotonic()
    def deadline():
        if time.monotonic()-tick>wall_seconds:raise TimeoutError('Interior source-sheet proof deadline')
    ids=np.asarray(face_ids,int)
    if ids.ndim!=1 or not len(ids) or len(np.unique(ids))!=len(ids) or ids.min()<0 or ids.max()>=len(mesh.faces):
        raise ValueError('Explicit distinct existing source facets required')
    v,alias=np.unique(mesh.vertices,axis=0,return_inverse=True);f=alias[mesh.faces]
    keep=np.setdiff1d(np.arange(len(f)),ids);cf=f[keep]
    olde,oldc=edges(f);newe,newc=edges(cf)
    before={tuple(e) for e in olde[oldc==1]};after={tuple(e) for e in newe[newc==1]}
    if not after<=before:raise ValueError('Source-sheet removal would open a new exterior boundary')
    local=topology(mesh.vertices,mesh.faces[ids])
    if local['edgeConnectedComponents']!=1 or local['boundaryLoops']!=1 or local['eulerCharacteristic']!=1:
        raise ValueError('The declared source patch must be one simple disk')
    candidate=Mesh(mesh.vertices,mesh.faces[keep],mesh.frame,mesh.attributes,mesh.source)
    d=diagnostics(candidate,include_pairs=True);deadline()
    if d['nonmanifold_edges'] or d['nonmanifold_vertices'] or d['self_intersection_pairs'] or not d['orientable']:
        raise ValueError('Retained complete source envelope must be manifold and self-disjoint')
    rings=boundary_loops(cf);caps=[];balls=[]
    for ring in rings:
        deadline();ring=np.asarray(ring,int)[::-1];points=v[ring];center=points.mean(0)
        _,_,basis=np.linalg.svd(points-center,full_matrices=False);xy=(points-center)@basis[:2].T
        fan=np.asarray(earcut(ring,{int(i):p for i,p in zip(ring,xy)}),int)
        if fan.shape!=(len(ring)-2,3):raise ValueError('Computational mouth closure lost a source control')
        caps.append(fan);balls.append((center,float(np.linalg.norm(points-center,axis=1).max())))
    cap_faces=np.concatenate(caps) if caps else np.empty((0,3),int)
    closed=Mesh(v,np.r_[cf,cap_faces],mesh.frame);dc=diagnostics(closed,include_pairs=True);deadline()
    if dc['nonmanifold_edges'] or dc['nonmanifold_vertices'] or dc['self_intersection_pairs'] or not dc['orientable'] or boundary_loops(closed.faces):
        raise ValueError('Separate computational closure is not a regular closed envelope')
    tri=mesh.triangles[ids];normal=np.cross(tri[:,1]-tri[:,0],tri[:,2]-tri[:,0]);area=np.linalg.norm(normal,axis=1)
    if np.any(area<=1e-20):raise ValueError('Degenerate source-sheet face')
    normal/=area[:,None];center=tri.mean(1);dist=BVH(candidate).closest_points(center)[0]['distance_m']
    if np.any(dist<=1e-7):raise ValueError('Declared sheet lies on the retained source exterior')
    epsilon=np.minimum(dist*.1,1e-5);probes=np.r_[center,center+normal*epsilon[:,None],center-normal*epsilon[:,None]]
    margin=min((float((np.linalg.norm(probes-c,axis=1)-r).min()) for c,r in balls),default=float('inf'))
    if margin<=1e-7:raise ValueError('Interior classification depends on an artificial mouth closure')
    origin=v.mean(0);w=igl.winding_number(v-origin,closed.faces,probes-origin);deadline()
    if np.max(np.abs(np.abs(w)-1))>1e-9:raise ValueError('Declared source sheet is not strictly inside the retained envelope')
    deleted=Mesh(mesh.vertices,mesh.faces[ids],mesh.frame)
    for other in [candidate,Mesh(v,cap_faces,mesh.frame)] if len(cap_faces) else [candidate]:
        _,contact=pair_contacts(deleted,other,0.,max_tests=200000);deadline()
        if not contact['complete_within_pair_budget'] or contact['crossing_face_pairs'] or contact['coplanar_overlap_face_pairs']:
            raise ValueError('Source partition intersects the retained exterior or a computational mouth cap')
    return keep,dict(complete=True,removedSourceFaceIds=ids.tolist(),retainedSourceFaceIds=keep.tolist(),
        sourcePatchTopology=local,retainedTopology=topology(mesh.vertices,mesh.faces[keep]),
        sourceOpenings=len(rings),removedFreeEdges=[list(x) for x in sorted(before-after)],
        newBoundaryEdges=0,allPositionsAndSurvivingSourceFacesLiteral=True,
        computationalCapFaces=cap_faces.tolist(),computationalCapsEmitted=0,
        minimumRemovedCentroidDistanceToExterior_m=float(dist.min()),
        minimumProbeClearanceFromMouthHullBalls_m=margin if np.isfinite(margin) else None,
        maximumInteriorWindingError=float(np.abs(np.abs(w)-1).max()),
        anatomicalArtifactIdentityInferred=False,anatomicalAcceptance=False)
