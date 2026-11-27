"""Remove an explicit isolated triangular fin from an outer vessel envelope.

The retained wall must already have opposed faces along both attachment edges.
The third edge must be an isolated free edge, never part of a real mouth.
No positions or surviving faces change; unrelated defects remain separate.
"""
import numpy as np

def remove_source_open_fin(mesh,source_face_id):
    face_id=int(source_face_id)
    if not 0<=face_id<len(mesh.faces):raise ValueError('Explicit existing source triangle required')
    v,alias=np.unique(mesh.vertices,axis=0,return_inverse=True);f=alias[mesh.faces];tri=f[face_id]
    if len(set(tri))!=3:raise ValueError('Degenerate fin cannot be classified')
    all_edges=np.sort(f[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
    edges,inv,count=np.unique(all_edges,axis=0,return_inverse=True,return_counts=True);lookup={tuple(e):i for i,e in enumerate(edges)}
    boundary=[tuple(sorted(e)) for e in zip(tri,np.roll(tri,-1))]
    if sorted(int(count[lookup[e]]) for e in boundary)!=[1,3,3]:raise ValueError('Fin must have two triple attachments and one free edge')
    attached=[e for e in boundary if count[lookup[e]]==3];free=next(e for e in boundary if count[lookup[e]]==1)
    old_open={tuple(e) for e in edges[count==1]}
    if any(set(free)&set(e) for e in old_open if e!=free):raise ValueError('The free edge belongs to an external source-mouth boundary')
    neighbours=[]
    for e in attached:
        ids=np.setdiff1d(np.unique(np.flatnonzero(inv==lookup[e])//3),[face_id]);assert len(ids)==2
        directions=[1 if e in set(zip(t,np.roll(t,-1))) else -1 for t in f[ids]]
        if sorted(directions)!=[-1,1]:raise ValueError('Retained carrier lacks consistently opposed wall faces')
        neighbours.append(dict(physicalEdge=list(e),retainedSourceFaceIds=ids.tolist()))
    keep=np.setdiff1d(np.arange(len(f)),[face_id]);retained=f[keep]
    if not set(tri)<=set(np.unique(retained)):raise ValueError('A free branch vertex would be deleted')
    ae,ac=np.unique(np.sort(retained[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1),axis=0,return_counts=True)
    new_open={tuple(e) for e in ae[ac==1]}
    if new_open!=old_open-{free}:raise ValueError('An actual carrier opening changed')
    return keep,dict(complete=True,removedSourceFaceIds=[face_id],retainedSourceFaceIds=keep.tolist(),
        attachmentEdges=neighbours,isolatedFreeEdgePoints_m=v[list(free)].tolist(),
        noNewBoundaryEdges=True,allOtherSourceBoundaryEdgesLiteral=True,allSurvivingFacesAndFieldsLiteral=True,
        newFaces=0,newCaps=0,anatomicalStructureInferred=False,anatomicalAcceptance=False)
