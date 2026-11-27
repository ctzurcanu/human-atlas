"""Topology proof for an explicitly declared, nonanatomical source fin.

The caller must establish that the source face is an artifact, not a valve or
other tissue. This combinatorial test alone cannot make that anatomical choice.
No geometric tolerance or coordinate welding modifies the input mesh.
"""
from collections import defaultdict
import numpy as np


def remove_declared_fin(mesh,face_id):
    face_id=int(face_id)
    if not 0<=face_id<len(mesh.faces):raise ValueError('Existing explicitly declared source face required')
    _,physical=np.unique(mesh.vertices,axis=0,return_inverse=True);f=physical[mesh.faces]
    incidence=defaultdict(list)
    for i,face in enumerate(f):
        for a,b in zip(face,np.roll(face,-1)):incidence[tuple(sorted([int(a),int(b)]))].append(i)
    edges=[tuple(sorted([int(a),int(b)])) for a,b in zip(f[face_id],np.roll(f[face_id],-1))]
    valences=[len(incidence[e]) for e in edges]
    if len(set(f[face_id]))!=3 or sorted(valences)!=[1,3,3]:
        raise ValueError('Source face is not an isolated free-edge fin attached along exactly two triple edges')
    for e in edges:
        survivors=[i for i in incidence[e] if i!=face_id]
        if len(survivors)==2:
            directions=[]
            for i in survivors:
                cyclic=list(zip(f[i],np.roll(f[i],-1)))
                directions.append(1 if e in cyclic else -1)
            if directions[0]==directions[1]:raise ValueError('Fin removal would leave inconsistent surface winding')
    keep=np.delete(np.arange(len(f)),face_id)
    before={e for e,ids in incidence.items() if len(ids)==1}
    after={e for e,ids in incidence.items() if len(ids)-(face_id in ids)==1}
    if not after<=before:raise ValueError('Fin removal must not create a new boundary edge')
    return keep,dict(complete=True,declaredSourceFace=face_id,edgeValencesBefore=valences,
        fixedSourceCoordinates=True,newBoundaryEdges=[],removedFreeBoundaryEdges=[list(e) for e in sorted(before-after)],
        tripleEdgesBecomeConsistentlyOrientedTwoFaceEdges=True,retainedSourceFaceIds=keep.tolist(),
        anatomicalArtifactIdentityInferred=False,anatomicalAcceptance=False)
