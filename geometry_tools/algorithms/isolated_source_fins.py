"""Identify an isolated two-triangle sheet attached to one regular carrier edge.

This is deliberately not a size-based mesh cleanup. A candidate must be an open
quad sharing exactly one edge and its two vertices with a larger source sheet.
The larger sheet already has both opposed incident faces along that edge. No
tube wall is filled, displaced, retriangulated or replaced.
"""
import numpy as np
from scipy import sparse
from scipy.sparse.csgraph import connected_components


def isolated_source_fins(mesh):
    v,alias=np.unique(mesh.vertices,axis=0,return_inverse=True)
    f=alias[mesh.faces]
    edges=np.sort(f[:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
    unique,index,count=np.unique(edges,axis=0,return_inverse=True,return_counts=True)
    rr=[];cc=[]
    for ei in np.flatnonzero(count==2):
        faces=np.flatnonzero(index==ei)//3
        rr.extend(faces.tolist());cc.extend(faces[::-1].tolist())
    _,labels=connected_components(sparse.coo_matrix((np.ones(len(rr)),(rr,cc)),shape=(len(f),len(f))).tocsr())
    patches=[];removed=set()
    for label in np.unique(labels):
        ids=np.flatnonzero(labels==label)
        if len(ids)!=2:continue
        vertices=np.unique(f[ids]);local=np.sort(f[ids][:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
        le,lc=np.unique(local,axis=0,return_counts=True)
        if len(vertices)!=4 or len(le)!=5 or sorted(lc.tolist())!=[1,1,1,1,2]:continue
        lookup={tuple(e):int(k) for k,e in enumerate(unique)}
        boundary=[lookup[tuple(e)] for e in le[lc==1]]
        attached=[k for k in boundary if count[k]!=1]
        if len(attached)!=1 or count[attached[0]]!=3:continue
        edge=unique[attached[0]]
        neighbours=np.setdiff1d(np.flatnonzero(index==attached[0])//3,ids)
        if len(neighbours)!=2 or labels[neighbours[0]]!=labels[neighbours[1]]:continue
        carrier=np.flatnonzero(labels==labels[neighbours[0]])
        if len(carrier)<16:continue
        if set(vertices)&set(np.unique(f[carrier]))!=set(edge):continue
        # The retained wall is already consistently oriented across its edge.
        directions=[]
        for face in f[neighbours]:
            cycle=list(zip(face,np.roll(face,-1)))
            directions.append(1 if tuple(edge) in cycle else -1 if tuple(edge[::-1]) in cycle else 0)
        if sorted(directions)!=[-1,1]:continue
        patches.append(dict(sourceFaceIds=ids.tolist(),sharedPhysicalEdge=edge.tolist(),
            sharedEdgePoints_m=v[edge].tolist(),retainedIncidentSourceFaceIds=neighbours.tolist(),
            carrierSourceFaceCount=len(carrier),freeBoundaryEdges=3,removedSourceSheetIsAnOpenQuad=True,
            retainedCarrierAlreadyHasOpposedEdgeFaces=True))
        removed.update(map(int,ids))
    removed=np.array(sorted(removed),int);kept=np.setdiff1d(np.arange(len(f)),removed)
    # Removing a fin is not allowed to create even one new carrier boundary.
    after_edges=np.sort(f[kept][:,[[0,1],[1,2],[2,0]]].reshape(-1,2),axis=1)
    ae,ac=np.unique(after_edges,axis=0,return_counts=True)
    old_open={tuple(e) for e in unique[count==1]};new_open={tuple(e) for e in ae[ac==1]}
    assert new_open<=old_open
    return kept,dict(complete=True,removedSourceFaceIds=removed.tolist(),retainedSourceFaceIds=kept.tolist(),
        patches=patches,nonmanifoldEdgesBefore=int((count>2).sum()),nonmanifoldEdgesAfter=int((ac>2).sum()),
        boundaryEdgesBefore=int((count==1).sum()),boundaryEdgesAfter=int((ac==1).sum()),
        newBoundaryEdges=0,allRetainedTrianglesLiteral=True,positionsAndAttributesModified=False,
        anatomicalStructureInferred=False)
