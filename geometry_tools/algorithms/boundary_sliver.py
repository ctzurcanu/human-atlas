"""Explicit microscopic boundary-sliver contraction, never a hole filler.

The caller names a real source face. A unique short edge must connect a
degree-two sliver point to the offending three-face edge. All other controls
and remaining source faces are literal. Anatomical acceptance is external.
"""
import numpy as np
from ..model import Mesh
from .tube_continuity import boundary_loops


def collapse_boundary_sliver(mesh,face_id,*,maximum_motion_m):
    if not isinstance(face_id,int)or not 0<=face_id<len(mesh.faces)or not 0<maximum_motion_m<=.00001:
        raise ValueError('Explicit source face and contraction bound up to10micrometres required')
    v,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inverse[mesh.faces]
    edges=np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1)
    unique,ids,count=np.unique(edges,axis=0,return_inverse=True,return_counts=True)
    face=f[face_id];face_edges=np.sort(np.c_[face,np.roll(face,-1)],axis=1)
    bad=[edge for edge in face_edges if count[np.flatnonzero(np.all(unique==edge,axis=1))[0]]==3]
    if len(bad)!=1:raise ValueError('The selected sliver must touch exactly one three-face source edge')
    lengths=np.linalg.norm(v[face]-v[np.roll(face,-1)],axis=1);short=np.flatnonzero(lengths<=maximum_motion_m)
    if len(short)!=1:raise ValueError('One unique microscopic source edge required')
    edge=face_edges[short[0]];ends=set(map(int,bad[0]));inside=[int(i)for i in edge if i in ends];outside=[int(i)for i in edge if i not in ends]
    if len(inside)!=1 or len(outside)!=1:raise ValueError('Short edge must terminate on the offending source edge')
    target,source=inside[0],outside[0]
    incident=np.flatnonzero(np.any(f==source,axis=1))
    if len(incident)!=2:raise ValueError('The sliver tip must have exactly two incident source faces')
    moved=np.flatnonzero(inverse==source);q=mesh.vertices.copy();q[moved]=v[target]
    triangles=q[mesh.faces];area=np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1)
    zero=np.flatnonzero(area==0)
    if zero.tolist()!=[face_id]:raise ValueError('Contraction must collapse only the declared source sliver')
    retained=np.flatnonzero(area>0);new_faces=mesh.faces[retained]
    # All new boundary cycles must be simple and consistently wound. This
    # establishes source topology only, never a physiological open lumen.
    physical,inv=np.unique(q,axis=0,return_inverse=True)
    loops=boundary_loops(inv[new_faces])
    before=mesh.triangles[retained];after=triangles[retained]
    normals=np.cross(before[:,1]-before[:,0],before[:,2]-before[:,0]);changed=np.cross(after[:,1]-after[:,0],after[:,2]-after[:,0])
    if np.any(np.einsum('ij,ij->i',normals,changed)<=0):raise ValueError('A retained source face inverted')
    displacement=float(np.linalg.norm(q-mesh.vertices,axis=1).max())
    return Mesh(q,new_faces,mesh.frame,mesh.attributes,mesh.source),dict(sourceFace=face_id,
        retainedSourceFaceIds=retained.tolist(),movedRawControls=moved.tolist(),sourcePhysicalPoint=source,
        fixedDestinationPhysicalPoint=target,sourcePoint_m=v[source].tolist(),destinationPoint_m=v[target].tolist(),
        maximumMotion_m=displacement,maximumMotionBound_m=maximum_motion_m,
        removedSourceArea_m2=float(np.linalg.norm(np.cross(mesh.triangles[face_id,1]-mesh.triangles[face_id,0],mesh.triangles[face_id,2]-mesh.triangles[face_id,0]))*.5),
        boundaryLoops=len(loops),boundaryLoopLengths=list(map(len,loops)),
        sourceDetailPreservedExceptDeclaredMicroscopicSliver=True,newFacesAdded=0,newCapsAdded=0,anatomicalAcceptance=False)


def remove_dangling_sliver_sheet(mesh,face_ids,*,maximum_width_m):
    """Remove an explicitly selected microscopic sheet attached by one edge.

    Nothing is welded or moved. The retained carrier must have simple source
    mouths and no nonmanifold edge. A normal vessel branch or broad surface
    cannot satisfy the one-edge attachment and microscopic-width guards.
    """
    selected=np.asarray(face_ids,int)
    if (selected.ndim!=1 or not 1<=len(selected)<=16 or len(np.unique(selected))!=len(selected)
        or np.any(selected<0)or np.any(selected>=len(mesh.faces))or not 0<maximum_width_m<=.00001):
        raise ValueError('Explicit bounded source-sheet faces and width required')
    v,inverse=np.unique(mesh.vertices,axis=0,return_inverse=True);f=inverse[mesh.faces]
    keep=np.flatnonzero(~np.isin(np.arange(len(f)),selected));strip=f[selected];retained=f[keep]
    shared=np.intersect1d(np.unique(strip),np.unique(retained))
    if len(shared)!=2:raise ValueError('The selected sheet must meet the retained carrier only at one edge')
    edge=np.sort(shared);all_edges=np.sort(f[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1)
    unique,count=np.unique(all_edges,axis=0,return_counts=True)
    hit=np.flatnonzero(np.all(unique==edge,axis=1))
    if len(hit)!=1 or count[hit[0]]!=3:raise ValueError('The sheet must be the third incident face at the retained vessel edge')
    strip_edges=np.sort(strip[:,[0,1,1,2,2,0]].reshape(-1,2),axis=1)
    e,ec=np.unique(strip_edges,axis=0,return_counts=True)
    if len(np.unique(strip))-len(e)+len(strip)!=1 or np.any(ec>2):raise ValueError('Selected source sheet is not a disk')
    graph={i:set()for i in range(len(strip))}
    for i in graph:
        for j in range(i):
            if len(set(strip[i])&set(strip[j]))==2:graph[i].add(j);graph[j].add(i)
    seen={0};todo=[0]
    while todo:
        for j in graph[todo.pop()]-seen:seen.add(j);todo.append(j)
    if len(seen)!=len(strip):raise ValueError('Selected source sheet is disconnected')
    exclusive=np.setdiff1d(np.unique(strip),shared)
    if not len(exclusive):raise ValueError('A retained vessel facet cannot be removed as a sliver sheet')
    distance=np.linalg.norm(v[exclusive,None]-v[shared][None],axis=2).min(1)
    if float(distance.max())>maximum_width_m:raise ValueError('The source sheet contains material beyond the microscopic endpoint width')
    loops=boundary_loops(retained)
    triangles=mesh.triangles[selected];area=float(np.linalg.norm(np.cross(triangles[:,1]-triangles[:,0],triangles[:,2]-triangles[:,0]),axis=1).sum()*.5)
    return Mesh(mesh.vertices.copy(),mesh.faces[keep].copy(),mesh.frame,
        {k:x.copy()for k,x in mesh.attributes.items()},mesh.source),dict(
        removedSourceFaceIds=selected.tolist(),retainedSourceFaceIds=keep.tolist(),
        attachmentPhysicalEdge=shared.tolist(),attachmentPoints_m=v[shared].tolist(),exclusivePhysicalPoints=exclusive.tolist(),
        maximumEndpointWidth_m=float(distance.max()),maximumWidthBound_m=maximum_width_m,removedSourceArea_m2=area,
        boundaryLoops=len(loops),boundaryLoopLengths=list(map(len,loops)),allSourceControlsLiteral=True,
        allSourceAttributesLiteral=True,newFacesAdded=0,newCapsAdded=0,anatomicalAcceptance=False)
