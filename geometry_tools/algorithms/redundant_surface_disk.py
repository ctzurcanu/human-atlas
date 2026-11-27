"""Remove an explicitly classified extra disk from an outer source envelope.

The disk must attach along a whole triple-incidence cycle. The retained source
already supplies both consistently wound faces on every cycle edge. A real
opening, annular branch, free valve leaflet or incomplete chart cannot qualify.
The caller must establish the outer-envelope representation and inspect it.
"""
import numpy as np
from geometry_tools.algorithms.shared_boundary_arc import topology
from geometry_tools.algorithms.tube_continuity import boundary_loops

def remove_redundant_surface_disk(mesh,source_face_ids):
    ids=np.asarray(source_face_ids,int)
    if ids.ndim!=1 or not len(ids) or len(np.unique(ids))!=len(ids) or ids.min()<0 or ids.max()>=len(mesh.faces):
        raise ValueError('Explicit distinct source disk faces required')
    v,alias=np.unique(mesh.vertices,axis=0,return_inverse=True);f=alias[mesh.faces]
    patch=f[ids];t=topology(mesh.vertices,mesh.faces[ids])
    if t['edgeConnectedComponents']!=1 or t['eulerCharacteristic']!=1 or t['boundaryLoops']!=1 or t['nonmanifoldEdges']:
        raise ValueError('Declared source patch is not one complete simple disk')
    ring=boundary_loops(patch)[0];boundary={tuple(sorted(e)) for e in zip(ring,np.roll(ring,-1))}
    all_edges=np.sort(f[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
    edges,inv,count=np.unique(all_edges,axis=0,return_inverse=True,return_counts=True);lookup={tuple(e):i for i,e in enumerate(edges)}
    if any(count[lookup[e]]!=3 for e in boundary):raise ValueError('A genuine opening or non-triple edge cannot be removed')
    keep=np.setdiff1d(np.arange(len(f)),ids);retained=f[keep]
    if set(np.unique(patch))&set(np.unique(retained))!=set(ring):raise ValueError('Disk interior shares another source structure')
    pairs=[]
    for e in sorted(boundary):
        owners=np.unique(np.flatnonzero(inv==lookup[e])//3);outside=np.setdiff1d(owners,ids)
        if len(outside)!=2:raise ValueError('Retained envelope does not already own two faces per attachment edge')
        directions=[]
        for face in f[outside]:
            cycle=set(zip(face,np.roll(face,-1)));directions.append(1 if e in cycle else -1)
        if sorted(directions)!=[-1,1]:raise ValueError('Removal would retain inconsistent source winding')
        pairs.append(dict(physicalEdge=list(e),retainedSourceFaces=outside.tolist()))
    after_edges,after_count=np.unique(np.sort(retained[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1),axis=0,return_counts=True)
    old_open={tuple(e) for e in edges[count==1]};new_open={tuple(e) for e in after_edges[after_count==1]}
    if old_open!=new_open:raise ValueError('A genuine source opening changed')
    return keep,dict(complete=True,removedSourceFaceIds=ids.tolist(),retainedSourceFaceIds=keep.tolist(),
        diskTopology=t,attachmentRingPhysicalPoints_m=v[ring].tolist(),attachmentEdges=pairs,
        sourceOpenBoundaryEdgesLiteral=True,allRetainedSourceTrianglesLiteral=True,
        sourcePositionsAndWordArraysLiteral=True,facesAdded=0,anatomicalIdentityInferred=False,anatomicalAcceptance=False)
