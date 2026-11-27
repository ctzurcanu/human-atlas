"""Recover explicit retained directed faces from a coordinate-replayed bank.

Geometry proposes original face identities. The final injective raw-control
map independently verifies them before source triangle pairs can be reused.
No missing original face is restored and no new quad partner is guessed.
"""
from collections import defaultdict,deque
import numpy as np
from scipy.spatial import cKDTree
from .source_face_pairing import remap_source_pairs


def retained_source_face_ancestry(source_vertices,source_faces,current_vertices,current_faces,*,maximum_error_m=1e-7):
    a,f,b,g=[np.asarray(x)for x in [source_vertices,source_faces,current_vertices,current_faces]]
    if not 0<maximum_error_m<=1e-7:raise ValueError('Explicit source replay precision within0.1micrometre required')
    physical,original_ids=np.unique(a,axis=0,return_inverse=True)
    distance,new_ids=cKDTree(physical).query(b)
    if float(distance.max(initial=0))>maximum_error_m:raise ValueError('Current controls do not replay the original source coordinate frame')
    lookup=defaultdict(deque)
    for i,face in enumerate(original_ids[f]):lookup[tuple(map(int,face))].append(i)
    kept=[]
    for face in new_ids[g]:
        key=tuple(map(int,face))
        if not lookup[key]:raise ValueError('A current directed face is not an unused original source face')
        kept.append(lookup[key].popleft())
    kept=np.asarray(kept,int)
    pairs,proof=remap_source_pairs(f,g,kept,np.arange(len(f)//2*2).reshape(-1,2))
    return kept,pairs,dict(proof,maximumCoordinateReplayError_m=float(distance.max(initial=0)),
                          removedOriginalFaceIds=np.setdiff1d(np.arange(len(f)),kept).tolist(),
                          currentFaceOrderRetained=True,sourceFrameReplayRequired=True)
