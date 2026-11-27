"""Certify a declared near-copy sleeve from shared source-root incidence.

The caller must establish alternative-source identity and material ownership.
No nearby anatomy is matched or welded, and no coordinate is modified.
"""
import numpy as np

def verify_sleeve_correspondence(mesh,copy_faces,retained_faces,root_raw_ids,position_codes,*,maximum_world_difference_m=1e-6,maximum_code_difference=1,maximum_states=10000):
    copy,retained=[np.asarray(x,int) for x in [copy_faces,retained_faces]];codes=np.asarray(position_codes)
    if not len(copy) or len(copy)!=len(retained) or codes.shape!=mesh.vertices.shape:raise ValueError('Equal source-face counts and original position-code words required')
    if np.intersect1d(copy,retained).size:raise ValueError('The two declared alternatives must be disjoint source primitives')
    _,representative,physical=np.unique(mesh.vertices,axis=0,return_index=True,return_inverse=True)
    first,second=[physical[mesh.faces[ids]] for ids in [copy,retained]];root=set(map(int,physical[np.asarray(root_raw_ids,int)]))
    if len(root)<3 or not root<=set(first.ravel()) or not root<=set(second.ravel()):raise ValueError('An explicit common source root is required')
    def graph(f):
        g={int(x):set() for x in np.unique(f)}
        for tri in f:
            for a,b in zip(tri,np.roll(tri,-1)):g[int(a)].add(int(b));g[int(b)].add(int(a))
        return g
    ga,gb=graph(first),graph(second)
    if len(ga)!=len(gb):raise ValueError('Source sleeve control topology differs')
    cycle=lambda t:tuple(map(int,np.roll(t,-int(np.argmin(t)))))
    target={cycle(t) for t in second};states=0
    def search(mapping):
        nonlocal states
        states+=1
        if states>maximum_states:raise ValueError('Finite rooted correspondence budget exceeded')
        if len(mapping)==len(ga):return mapping if {cycle([mapping[int(x)] for x in tri]) for tri in first}==target else None
        node=max(set(ga)-set(mapping),key=lambda n:len(ga[n]&mapping.keys()))
        for dest in sorted(set(gb)-set(mapping.values())):
            if len(ga[node])!=len(gb[dest]) or any((old in ga[node])!=(new in gb[dest]) for old,new in mapping.items()):continue
            result=search(mapping|{node:dest})
            if result is not None:return result
        return None
    mapping=search({x:x for x in root})
    if mapping is None:raise ValueError('No rooted oriented triangle correspondence exists')
    first_ids=np.array([representative[x] for x in mapping]);second_ids=np.array([representative[mapping[x]] for x in mapping]);delta=float(np.linalg.norm(mesh.vertices[first_ids]-mesh.vertices[second_ids],axis=1).max())
    words=np.abs(codes[first_ids].astype(np.int64)-codes[second_ids].astype(np.int64))
    if delta>maximum_world_difference_m or words.max()>maximum_code_difference:raise ValueError('Sleeve difference exceeds the declared original source-code/world domain')
    # Every oriented source triangle has a unique retained triangle with the
    # same barycentric parameter domain. The norm of the affine difference is
    # convex, hence its maximum is bounded by the maximum corner difference.
    return dict(complete=True,passed=True,physicalVertexCorrespondence={str(x):int(y) for x,y in mapping.items()},maximumMatchedWorldDifference_m=delta,maximumPositionCodeDifference=words.max(0).tolist(),symmetricSurfaceDistanceUpperBound_m=delta,orientedSourceFaceBijection=True,allSourceCornersIncluded=True,sourcePositionsModified=False,sourceFieldsModified=False,anatomicalIdentityInferred=False,correspondenceStates=states)
